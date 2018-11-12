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
    log = parser.getbool( step, 'log', False )

    parser.checkSpelling( step, soltab, ['axesToSmooth', 'size', 'mode', 'degree', 'replace', 'log'])
    return run(soltab, axesToSmooth, size, mode, degree, replace, log)

def _savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.

    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)

    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).

    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size -1 ) // 2

    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)

    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))

    # TODO: use astropy.convolution.convolve here to deal better with NaNs?
    return np.convolve( m[::-1], y, mode='valid')


def run( soltab, axesToSmooth, size=[], mode='runningmedian', degree=1, replace=False, log=False):
    """
    A smoothing function: running-median on an arbitrary number of axes, running polyfit and Savitzky-Golay on one axis, or set all solutions to the mean/median value.
    WEIGHT: flag ready.

    Parameters
    ----------
    axesToSmooth : array of str
        Axes used to compute the smoothing function.

    size : array of int, optional
        Window size for the runningmedian, savitzky-golay, and runningpoly (array of same size of axesToSmooth), by default [].

    mode : {'runningmedian','runningpoly','savitzky-golay','mean','median'}, optional
        Runningmedian or runningpoly or Savitzky-Golay or mean or median (these last two values set all the solutions to the mean/median), by default "runningmedian".

    degree : int, optional
        Degrees of the polynomia for the runningpoly or savitzky-golay modes, by default 1.

    replace : bool, optional
        Flagged data are replaced with smoothed value and unflagged, by default False.

    log : bool, optional
        clip is done in log10 space, by default False
    """

    import numpy as np
    from scipy.ndimage import generic_filter

    if mode == "runningmedian" and len(axesToSmooth) != len(size):
        logging.error("Axes and Size lengths must be equal for runningmedian.")
        return 1

    if (mode == "runningpoly" or mode=="savitzky-golay") and (len(axesToSmooth) != 1 or len(size) != 1):
        logging.error("Axes and size lengths must be 1 for runningpoly or savitzky-golay.")
        return 1

    if (mode == "runningpoly" or mode=="savitzky-golay") and soltab.getType() == 'phase':
        logging.error("Runningpoly and savitzky-golay modes cannot work on phases.")
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

    if soltab.getType() == 'amplitude' and not log:
        logging.warning('Amplitude solution tab detected and log=False. Amplitude solution tables should be treated in log space.')

    if mode == 'median' or mode == 'mean':
        vals = soltab.getValues(retAxesVals=False)
        if log: vals = np.log10(vals)
        weights = soltab.getValues(retAxesVals=False, weight=True)
        np.putmask(vals, weights==0, np.nan)
        idx_axes = [soltab.getAxesNames().index(axisToSmooth) for axisToSmooth in axesToSmooth]

        # handle phases by using a complex array
        if soltab.getType() == 'phase':
            vals = np.exp(1j*vals)

        if mode == 'median':
            vals[:] = np.nanmedian( vals, axis=idx_axes, keepdims=True)
        if mode == 'mean':
            logging.warning('Mean does not support NaN yet, use median if it is a problem.')
            vals[:] = np.mean( vals, axis=idx_axes, keepdims=True) # annoying np.nanmean does not accept axis=list!

        # go back to phases
        if soltab.getType() == 'phase':
            vals = np.angle(vals)

        # write back
        if log: vals = 10**vals
        soltab.setValues(vals)
        if replace:
            weights[ (weights == 0) ] = 1
            weights[ np.isnan(vals) ] = 0 # all the slice was flagged, cannot estrapolate value
            soltab.setValues(weights, weight=True)

    else:
        for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToSmooth, weight=True):

            # skip completely flagged selections
            if (weights == 0).all(): continue
            if log: vals = np.log10(vals)

            if mode == 'runningmedian':
                vals_bkp = vals[ weights == 0 ]

                # handle phases by using a complex array
                if soltab.getType() == 'phase':
                    vals = np.exp(1j*vals)

                    valsreal = np.real(vals)
                    valsimag = np.imag(vals)
                    np.putmask(valsreal, weights == 0, np.nan)
                    np.putmask(valsimag, weights == 0, np.nan)

                    # run generic_filter twice, once for real once for imaginary
                    valsrealnew = generic_filter(valsreal, np.nanmedian, size=size, mode='constant', cval=np.nan)
                    valsimagnew = generic_filter(valsimag, np.nanmedian, size=size, mode='constant', cval=np.nan)
                    valsnew = valsrealnew + 1j*valsimagnew # go back to complex
                    valsnew = np.angle(valsnew) # go back to phases

                else: # other than phases
                    np.putmask(vals, weights == 0, np.nan)
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
                    weights[ np.isnan(valsnew) ] = 0 # all the size was flagged cannot extrapolate value
                else:
                    valsnew[ weights == 0 ] = vals_bkp
                #print coord['ant'], vals, valsnew

            elif mode == 'savitzky-golay':
                vals_bkp = vals[ weights == 0 ]
                np.putmask(vals, weights==0, np.nan)
                valsnew = _savitzky_golay(vals, size[0], degree)
                if replace:
                    weights[ weights == 0] = 1
                    weights[ np.isnan(valsnew) ] = 0 # all the size was flagged cannot extrapolate value
                else:
                    valsnew[ weights == 0 ] = vals_bkp

            else:
                logging.error('Mode must be: runningmedian, runningpoly, savitzky-golay, median or mean')
                return 1

            if log: valsnew = 10**valsnew
            soltab.setValues(valsnew, selection)
            if replace: soltab.setValues(weights, selection, weight=True)

    soltab.flush()
    soltab.addHistory('SMOOTH (over %s with mode = %s)' % (axesToSmooth, mode))
    return 0


