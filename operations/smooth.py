#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a smoothing function
#

import numpy as np
import logging
import scipy.ndimage.filters
from operations_lib import *
from h5parm import solFetcher, solWriter

logging.debug('Loading SMOOTH module.')

def run( step, parset, H ):

    soltabs = getParSoltabs( step, parset, H )
    ants = getParAnts( step, parset, H )
    pols = getParPols( step, parset, H )
    dirs = getParDirs( step, parset, H )

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

        for vals, coord, nrows in sf.getIterValuesGrid(returnAxes=axesToSmooth, return_nrows=True):
            # TODO: implement flag control
            # smoothing
            valsnew = scipy.ndimage.filters.median_filter(vals, FWHM)

            # writing back the solutions
            print "write started"
            sw.setValuesGrid(valsnew, nrows)
            print "write ended"

#         sw.addHistory('SMOOTH (over %s with box size = %s, ants = %s, '
#             'pols = %s, dirs = %s)' % (axesToSmooth, FWHM, ants, pols, dirs))
        selection = sf.selection
        sw.addHistory('SMOOTH (over %s with box size = %s and selection = [%s])' % (axesToSmooth, FWHM, selection))
    return 0


