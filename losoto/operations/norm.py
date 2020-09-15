#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading NORM module.')

def _run_parser(soltab, parser, step):
    axesToNorm = parser.getarraystr( step, 'axesToNorm' ) # no default
    normVal = parser.getfloat( step, 'normVal', 1.)
    log = parser.getbool( step, 'log', False )

    parser.checkSpelling( step, soltab, ['axesToNorm','normVal','log'])
    return run(soltab, axesToNorm, normVal, log)

def run( soltab, axesToNorm, normVal = 1., log = False ):
    """
    Normalize the solutions to a given value
    WEIGHT: Weights compliant

    Parameters
    ----------
    axesToNorm : array of str
        Axes along which compute the normalization.

    normVal : float, optional
        Number to normalize to vals = vals * (normVal/valsMean), by default 1.

    log : bool, optional
        clip is done in log10 space, by default False.
    """
    import numpy as np
    
    logging.info("Normalizing soltab: "+soltab.name)

    # input check
    axesNames = soltab.getAxesNames()
    for normAxis in axesToNorm:
        if normAxis not in axesNames:
            logging.error('Normalization axis '+normAxis+' not found.')
            return 1

    if soltab.getType() == 'amplitude' and not log:
        logging.warning('Amplitude solution tab detected and log=False. Amplitude solution tables should be treated in log space.')

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToNorm, weight = True):

        if log: vals = np.log10(vals)

        # rescale solutions
        if np.all(weights == 0): continue # skip flagged selections
        valsMean = np.nanmean(vals[weights != 0])

        if log: vals[weights != 0] += np.log10(normVal)-valsMean
        else: vals[weights != 0] *= normVal/valsMean

        logging.debug("Rescaling by: "+str(normVal/valsMean))

        # writing back the solutions
        if log: vals = 10**vals
        soltab.setValues(vals, selection)

    soltab.flush()
    soltab.addHistory('NORM (on axis %s)' % (axesToNorm))

    return 0
