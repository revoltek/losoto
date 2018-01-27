#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

logging.debug('Loading NORM module.')

def _run_parser(soltab, parser, step):
    axesToNorm = parser.getarraystr( step, 'axesToNorm' ) # no default
    normVal = parser.getfloat( step, 'normVal', 1.)
    return run(soltab, axesToNorm, normVal)

def run( soltab, axesToNorm, normVal = 1. ):
    """
    Normalize the solutions to a given value
    WEIGHT: Weights compliant

    Parameters
    ----------
    axesToNorm : array of str
        Axes along which compute the normalization.

    normVal : float, optional
        Number to normalize to vals = vals * (normVal/valsMean), by default 1.
    """
    import numpy as np
    
    logging.info("Normalizing soltab: "+soltab.name)

    # input check
    axesNames = soltab.getAxesNames()
    for normAxis in axesToNorm:
        if normAxis not in axesNames:
            logging.error('Normalization axis '+normAxis+' not found.')
            return 1

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToNorm, weight = True):

        # rescale solutions
        if np.sum(weights) == 0: continue # skip flagged selections
        valsMean = np.average(vals, weights=weights)
        vals[weights != 0] *= normVal/valsMean
        logging.debug(str(coord))
        logging.debug("Rescaling by: "+str(normVal/valsMean))

        # writing back the solutions
        soltab.setValues(vals, selection)

    soltab.flush()
    soltab.addHistory('NORM (on axis %s)' % (axesToNorm))

    return 0
