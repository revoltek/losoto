#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a smoothing function for the clock
# WEIGHH: weight ready

import logging
from losoto.operations_lib import *

logging.debug('Loading SMOOTHCLOCK module.')

def run( step, parset, H ):

    import scipy.ndimage.filters
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    s = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "smooth"]), 10 )

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Smoothing soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        if sf.getType() != 'clock':
            logging.error('Only clock-type solutions can be run in this operation.')
            return 1

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        for vals, weights, coord, selection in sf.getValuesIter(returnAxes='time', weight=True):

            x=coord['time'][weights != 0]
            y=vals[weights != 0]
            weights = weights[weights != 0]
            spline = scipy.interpolate.UnivariateSpline(x, y, weights, k=1, s=1e-15)

            plot = True
            if plot:
                import matplotlib as mpl
                mpl.use("Agg")
                import matplotlib.pyplot as plt
                plt.plot(x, y, 'k.')
                plt.plot(x, spline(x), 'r-')
                plt.savefig('test.png')
                sys.exit(1)

            sw.selection = selection
            sw.setValues(spline(coord['time']))

        sw.addHistory('Smoothed with SMOOTHCLOCK.')
        del sf
        del sw
    return 0


