#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a smoothing function
# WEIGHH: not ready

import logging
from losoto.operations_lib import *

logging.debug('Loading SMOOTH module.')

def run( step, parset, H ):

    import scipy.ndimage.filters
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    axesToSmooth = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    FWHM = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "FWHM"]), [] )
    mode = parset.getString('.'.join(["LoSoTo.Steps", step, "Mode"]), "runningmedian" )

    if mode == "runningmedian" and len(axesToSmooth) != len(FWHM):
        logging.error("Axes and FWHM lenghts must be equal.")
        return 1

    if FWHM != [] and mode != "runningmedian":
        logging.warning("FWHM makes sense only with runningmedian mode, ignoring it.")

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Smoothing soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab, useCache = True) # remember to flush!

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        for i, axis in enumerate(axesToSmooth[:]):
            if axis not in sf.getAxesNames():
                del axesToSmooth[i]
                del FWHM[i]
                logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

        for vals, coord, selection in sf.getValuesIter(returnAxes=axesToSmooth):

            if mode == 'runningmedian':
                valsnew = scipy.ndimage.filters.median_filter(vals, FWHM)
            elif mode == 'median':
                valsnew = np.median(vals)
            elif mode == 'mean':
                valsnew = np.mean(vals)
            else:
                logging.error('Mode must be: runningmedian, median or mean')
                return 1

            sw.selection = selection
            sw.setValues(valsnew)

        sw.flush()
        sw.addHistory('SMOOTH (over %s with mode = %s)' % (axesToSmooth, mode))
    return 0


