#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

logging.debug('Loading SMOOTH module.')

def _run_parser(soltab, parser, step):
    axesToSmooth = parser.getarraystr( step, 'axesToSmooth' ) # no default
    size = parser.getarrayint( step, 'size', [] )
    mode = parser.getstr( step, 'mode', 'runningmedian' )
    degree = parser.getint( step, 'degree', 1 )
    replace = parser.getbool( step, 'replace', False )
    return run(soltab, axesToSmooth, size, mode, degree, replace)

def run( soltab, axesToSmooth, size=[], mode='runningmedian', degree=1, replace=False):
    """
    This operation for LoSoTo implement a smoothing function: running-median on an arbitrary number of axes, running polyfit on one axis, or set all solutions to mean/median value.
    WEIGHH: flag ready

    Parameters
    ----------
    axesToSmooth : array of str
        Axes used to compute the smoothing function.
    
    size : array of int, optional
        Window size for the runningmedian and runningpoly (array of same size of axesToSmooth), by default [].

    mode : {'runningmedian','runningpoly','mean','median'}, optional
        Runningmedian or runningpoly or mean or median (these last two values set all the solutions to the mean/median), by default "runningmedian".

    degree : int, optional
        Degrees of the polynomia for the runningpoly, by default 1.

    replace : bool, optional
        Flagged data are replaced with smoothed value and unflagged, by default False.
    """

    import numpy as np
    from scipy.ndimage import generic_filter

    if mode == "runningmedian" and len(axesToSmooth) != len(size):
        logging.error("Axes and Size lenghts must be equal for runningmedian.")
        return 1

    if mode == "runningpoly" and (len(axesToSmooth) != 1 or len(size) != 1):
        logging.error("Axes and size lenghts must be 1 for runningpoly.")
        return 1

    for i, s in enumerate(size):
        if s % 2 == 0:
            logging.warning('Size should be odd, adding 1.')
            size[i] += 1

    logging.info("Smoothing soltab: "+soltab.name)

    for i, axis in enumerate(axesToSmooth[:]):
        if axis not in soltab.getAxesNames():
            del axesToSmooth[i]
            del size[i]
            logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToSmooth, weight=True):

        # skip completely flagged selections
        if (weights == 0).all(): continue

        if mode == 'runningmedian':
            vals_bkp = vals[ weights == 0 ]
            np.putmask(vals, weights==0, np.nan)
            valsnew = generic_filter(vals, np.nanmedian, size=size, mode='constant', cval=np.nan)
            if replace: 
                weights[ weights == 0] = 1
                weights[ np.isnan(valsnew) ] = 0 # all the size was flagged cannoth estrapolate value
            else:
                valsnew[ weights == 0 ] = vals_bkp

        elif mode == 'runningpoly':
            def polyfit(data):
                if (np.isnan(data)).all(): return np.nan # all size is flagged
                x = np.arange(len(data))[ ~np.isnan(data)]
                y = data[ ~np.isnan(data) ]
                p = np.polynomial.polynomial.polyfit(x, y, deg=degree)
                #import matplotlib as mpl
                #mpl.use("Agg")
                #import matplotlib.pyplot as plt
                #plt.plot(x, y, 'ro')
                #plt.plot(x, np.polyval( p[::-1], x ), 'k-')
                #plt.savefig('test.png')
                #sys.exit()
                return np.polyval( p[::-1], (size[0]-1)/2 ) # polyval has opposite convention for polynomial order

            # flags and at edges pass 0 and then remove them
            vals_bkp = vals[ weights == 0 ]
            np.putmask(vals, weights==0, np.nan)
            valsnew = generic_filter(vals, polyfit, size=size[0], mode='constant', cval=np.nan)
            if replace:
                weights[ weights == 0] = 1
                weights[ np.isnan(valsnew) ] = 0 # all the size was flagged cannoth estrapolate value
            else:
                valsnew[ weights == 0 ] = vals_bkp
            #print coord['ant'], vals, valsnew

        # TODO: if phase use angular mean/median
        elif mode == 'median':
            valsnew = np.nanmedian( vals[(weights!=0)] )
            if replace: weights[ weights == 0] = 1

        elif mode == 'mean':
            valsnew = np.nanmean( vals[(weights!=0)] )
            if replace: weights[ weights == 0] = 1

        else:
            logging.error('Mode must be: runningmedian, runningpoly, median or mean')
            return 1

        soltab.setValues(valsnew, selection)
        if replace: soltab.setValues(weights, selection, weight=True)

    soltab.flush()
    soltab.addHistory('SMOOTH (over %s with mode = %s)' % (axesToSmooth, mode))
    return 0


