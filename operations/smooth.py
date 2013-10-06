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
    ants = getParAnts( step, parset, H )
    pols = getParPols( step, parset, H )
    dirs = getParDirs( step, parset, H )
    logging.debug("Time for getting parms: "+str(elapsed)+" s.")

    axesToSmooth = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    FWHM = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "FWHM"]), [] )

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Smoothing soltab: "+soltab.name)

        sf.makeSelection(ant=ants, pol=pols, dir=dirs)

        # some checks
        if len(FWHM) != len(axesToSmooth):
            logging.error('Wrong length of FWHM/Axes parameter.')
            return 1

        for i, axis in enumerate(axesToSmooth[:]):
            if axis not in sf.getAxes():
                del axesToSmooth[i]
                del FWHM[i]
                logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

        start = time.clock()
        first = True
        for vals, coord, nrows in sf.getIterValuesGrid(returnAxes=axesToSmooth, return_nrows=True):
            if first == True:
                elapsed = (time.clock() - start)
                logging.debug("Time for load: "+str(elapsed)+" s.")
                first = False

            # TODO: implement flag control
            # smoothing
            valsnew = scipy.ndimage.filters.median_filter(vals, FWHM)
        
            # writing back the solutions
            sw.setValuesGrid(valsnew, nrows)
        start = time.clock()
        sw.flush()
        elapsed = (time.clock() - start)
        logging.debug("Time for write: "+str(elapsed)+" s.")

#         sw.addHistory('SMOOTH (over %s with box size = %s, ants = %s, '
#             'pols = %s, dirs = %s)' % (axesToSmooth, FWHM, ants, pols, dirs))
        selection = sf.selection
        sw.addHistory('SMOOTH (over %s with box size = %s and selection = [%s])' % (axesToSmooth, FWHM, selection))
    return 0


