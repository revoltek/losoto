#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an interpolation script for LoSoTo

import logging
from operations_lib import *

logging.debug('Loading NORM module.')

def run( step, parset, H ):
    """
    Normalize the solutions to a given value
    """
    import numpy as np
    from h5parm import solFetcher, solWriter
    
    solsets = getParSolsets( step, parset, H )
    soltabs = getParSoltabs( step, parset, H )
    solTypes = getParSolTypes( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    normVal = parset.getFloat('.'.join(["LoSoTo.Steps", step, "NormVal"]), 1. )
    normAxis = parset.getString('.'.join(["LoSoTo.Steps", step, "NormAxis"]), 'time' )

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab._v_name)

        tr = solFetcher(soltab)
        tw = solWriter(soltab)

        axesNames = tr.getAxesNames()
        if normAxis not in axesNames:
            logging.error('Normalization axis '+normAxis+' not found.')
            return 1

        tr.setSelection(ant=ants, pol=pols, dir=dirs)
        for vals, coord in tr.getValuesIter(returnAxes=normAxis):

            # construct grid
            coordSel = removeKeys(coord, normAxis)
            logging.debug("Working on coords:"+str(coordSel))

            # rescale solutions
            valsMean = np.mean(vals)
            valsNew = normVal*vals/valsMean
            logging.debug("Rescaling by: "+str(normVal/valsMean))

            # writing back the solutions
            tw.setSelection(**coordSel)
            tw.setValues(valsNew)

    tw.addHistory('NORM (on axis %s)' % (normAxis))
    return 0
