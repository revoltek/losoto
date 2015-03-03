#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an interpolation script for LoSoTo
# WEIGHT: Weights compliant

import logging
from operations_lib import *

logging.debug('Loading NORM module.')

def run( step, parset, H ):
    """
    Normalize the solutions to a given value
    """
    import numpy as np
    from h5parm import solFetcher, solWriter
    
    soltabs = getParSoltabs( step, parset, H )
    solTypes = getParSolTypes( step, parset, H )

    normVal = parset.getFloat('.'.join(["LoSoTo.Steps", step, "NormVal"]), 1. )
    normAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "NormAxes"]), ['time'] )

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Normalizing soltab: "+soltab._v_name)

        tr = solFetcher(soltab)
        tw = solWriter(soltab, useCache = True) # remember to flush!

        axesNames = tr.getAxesNames()
        for normAxis in normAxes:
            if normAxis not in axesNames:
                logging.error('Normalization axis '+normAxis+' not found.')
                return 1

        # axis selection
        userSel = {}
        for axis in tr.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        tr.setSelection(**userSel)

        for vals, weights, coord, selection in tr.getValuesIter(returnAxes=normAxes, weight = True):

            # rescale solutions
            if np.sum(weights) == 0: continue # skip flagged antenna
            valsMean = np.average(vals, weights=weights)
            valsNew = normVal*vals/valsMean
            logging.debug(str(coord))
            logging.debug("Rescaling by: "+str(normVal/valsMean))

            # writing back the solutions
            tw.selection = selection
            tw.setValues(valsNew)

        tw.flush()
        tw.addHistory('NORM (on axis %s)' % (normAxes))

    return 0
