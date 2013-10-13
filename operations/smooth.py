#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a smoothing function
#

import logging
from operations_lib import *

logging.debug('Loading SMOOTH module.')

def run( step, parset, H ):

    import time 

    import numpy as np
    import scipy.ndimage.filters
    from h5parm import solFetcher, solWriter

    start = time.clock()
    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )
    elapsed = (time.clock() - start)
    logging.debug("Time for getting parms: "+str(elapsed)+" s.")

    axesToSmooth = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    FWHM = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "FWHM"]), [] )

    if len(axesToSmooth) != len(FWHM):
        logging.error("Axes and FWHM lenghts must be equal.")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Smoothing soltab: "+soltab.name)

        sf.setSelection(ant=ants, pol=pols, dir=dirs)

        for i, axis in enumerate(axesToSmooth[:]):
            if axis not in sf.getAxesNames():
                del axesToSmooth[i]
                del FWHM[i]
                logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

        start = time.clock()
        first = True
        for vals, coord in sf.getValuesIter(returnAxes=axesToSmooth):
            if first == True:
                elapsed = (time.clock() - start)
                logging.debug("Time for load: "+str(elapsed)+" s.")
                first = False

            # TODO: implement flag control
            valsnew = scipy.ndimage.filters.median_filter(vals, FWHM)
        
            # writing back the solutions
            sw.setSelection(*coord) # TODO: do not select on returnAxes!!! much faster
            sw.setValuesGrid(valsnew)

        sw.addHistory('SMOOTH (over %s with box size = %s)' % (axesToSmooth, FWHM))
    return 0


