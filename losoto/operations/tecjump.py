#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a jumps remover for TEC solutions
# WEIGHH: flag ready

import logging
from losoto.operations_lib import *

logging.debug('Loading TECJUMP module.')

def run( step, parset, H ):

    import scipy.ndimage.filters
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    def robust_std(data, sigma=3):
        """
        Calculate standard deviation excluding outliers
        ok with masked arrays
        """
        return np.std(data[np.where(np.abs(data) < sigma * np.std(data))])

    def mask_interp(vals, mask):
        """
        return interpolated values for masked elements
        """
        this_vals = vals.copy()
        this_vals[mask] = np.interp(np.where(mask)[0], np.where(~mask)[0], vals[~mask])
        return this_vals

    def rolling_std(a, window, robust=False):
        """
        Return the rms for each element of the array calculated using the element inside a window,
        edges are mirrored
        """
        assert window % 2 == 1 # window must be odd
        a = np.pad(a, window/2, mode='reflect')
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        return np.sqrt(np.var(np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides), -1))

    win = 11
    nsigma = 5

    soltabs = getParSoltabs( step, parset, H )

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Smoothing soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab) # remember to flush!

        # TODO: check if it's a Tec table

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        for vals, weights, coord, selection in sf.getValuesIter(returnAxes='time', weight=True):

            # skip all flag
            if (weights == 0).all(): continue

            # kill large values
            weights[vals>0.5] = 0
            # interpolate flagged values to get resonable distances
            vals = mask_interp(vals, mask=(weights == 0))
            # get tev[i] - tec[i+1] - len is len(vals)-1
            vals_d = vals[:-1] - vals[1:]
            # skip reference
            if (vals_d == 0).all(): continue

            # get rolling std of distances
            std_d = rolling_std( mask_interp(vals_d, mask=(abs(vals_d)>0.01)), 51 )
            # get smooth distances
            smooth_d = scipy.ndimage.filters.median_filter( mask_interp(vals_d, mask=(abs(vals_d)>0.01)), 11 )
    
            f_dist = []
            jumps_init = []
            idx_jumps = []
            for i, d in enumerate(vals_d):
                if np.abs(d) > 5*std_d[i]:
                    idx_jumps.append(i)
                    jumps_init.append(d - smooth_d[i])

            #print vals_d[np.where(vals_d!=0)]
            print "%s: number of jumps: %i" % (coord['ant'], len(np.where(vals_d != 0)[0]))

            # couple of idexes for contiguos regions
            idx_jumps = zip([0]+idx_jumps,idx_jump+[len(vals_d)])

            for i, j in idx_jumps:
                    f_dist.append(lambda)

            # minimise the distance between each point and the std_d having vals_d as free parameters
            # define system of equation
            def f_all(f_dist, x)
                return np.sum(f(x[i]) for i,f in enumerate(f_dist))

            # correct vals with cumulative jumps
            #for i in range(len(vals)):
            #    vals[i] += np.sum(vals_d[0:i])

            # set back to 0 the values for flagged data
            vals[weights == 0] = 0

            sw.selection = selection
            sw.setValues(vals)
            sw.setValues(weights, weight=True)

        sw.addHistory('TECJUMP')
        del sf
        del sw
    return 0


