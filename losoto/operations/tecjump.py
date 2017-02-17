#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a jumps remover for TEC solutions
# WEIGHH: flag ready

import logging
from losoto.operations_lib import *

logging.debug('Loading SMOOTH module.')

def run( step, parset, H ):

    import scipy.ndimage.filters
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    def robust_std(data, sigma=3):
        """
        Calculate standard deviation excluding outliers
        ok with masked arrays
        """
        return np.std(data_d[np.where(np.abs(data_d) < sigma * np.std(data_d))])

    def mask_interp(vals, mask):
        """
        return interpolated values for masked elements
        """
        vals[mask] = np.interp(np.where(mask)[0], np.where(~mask)[0], vals[~mask])
        return vals

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

    win = 21
    nsigma = 3

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

        for vals, weights, coord, selection in sf.getValuesIter(returnAxes='ant', weight=True):

            # interpolate flagged values to get resonable distances
            vals = mask_interp(vals, mask=(weight == 0))
            # get tev[i] - tec[i+1] - len is len(vals)-1
            vals_d = vals[:-2] - vals[1:]
            # get rolling std of distances - TODO: replace high std vals_d with interpolated values?
            std_d = rolling_rms(vals_d, win)
            # get smooth distances - TODO: replace high std vals_d with interpolated values?
            smooth_d = scipy.ndimage.filters.median_filter(vals_d, FWHM)

            for i, d in enumerate(vals_d):
                if d < nsigma*std_d[i]: vals_d = 0 # no jump, leave 0
                # jump, replace distance with current distance minus expected distance (from smooth),
                # this should estimate the jump value
                vals_d[i] = d - smooth_d[i]

            # correct vals with cumulative jumps
            for i in range(len(vals)):
                vals[i] += np.sum(vals_d[0:i])


            # set back to 0 the values for flagged data
            vals[weight == 0] = 0

            sw.selection = selection
            sw.setValues(vals)

        sw.addHistory('SMOOTH (over %s with mode = %s)' % (axesToSmooth, mode))
        del sf
        del sw
    return 0


