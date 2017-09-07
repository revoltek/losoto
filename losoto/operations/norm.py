#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.operations_lib import *

logging.debug('Loading NORM module.')

def run_parser(soltab, parser, step):
    normAxes = parser.getarray( step, 'normAxes' ) # no default
    normVal = parser.getfloat( step, 'normVal', 1.)
    return run(soltab, normAxes, normVal)

def run( soltab, normAxes, normVal = 1. ):
    """
    Normalize the solutions to a given value
    WEIGHT: Weights compliant

    Parameters
    ----------
    normAxes : array of str
        Axes along which compute the normalization.

    normVal : float, optional
        Number to normalize to vals = vals * (normVal/valsMean), by default 1.
    """
    import numpy as np
    
    logging.info("Normalizing soltab: "+soltab.name)

    # input check
    axesNames = soltab.getAxesNames()
    for normAxis in normAxes:
        if normAxis not in axesNames:
            logging.error('Normalization axis '+normAxis+' not found.')
            return 1

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=normAxes, weight = True):

        # rescale solutions
        if np.sum(weights) == 0: continue # skip flagged selections
        valsMean = np.average(vals, weights=weights)
        vals[weights != 0] *= normVal/valsMean
        logging.debug(str(coord))
        logging.debug("Rescaling by: "+str(normVal/valsMean))

        # writing back the solutions
        soltab.setValues(vals, selection)

    soltab.flush()
    soltab.addHistory('NORM (on axis %s)' % (normAxes))

    return 0
