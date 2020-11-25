#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading SMOOTH module.')

def _run_parser(soltab, parser, step):
    axesToSmooth = parser.getarraystr( step, 'axesToSmooth' ) # no default
    size = parser.getarrayint( step, 'size', [] )
    mode = parser.getstr( step, 'mode', 'runningmedian' )
    degree = parser.getint( step, 'degree', 1 )
    replace = parser.getbool( step, 'replace', False )
    log = parser.getbool( step, 'log', False )
    refAnt = parser.getstr( step, 'refAnt', '' )

    parser.checkSpelling( step, soltab, ['axesToSmooth', 'size', 'mode', 'degree', 'replace', 'log', 'refAnt'])
    return run(soltab, axesToSmooth, size, mode, degree, replace, log, refAnt)


def _savitzky_golay(y, window_size, order):
    from scipy.signal import savgol_filter

    # replace any NaNs using linear interpolation
    nans = np.isnan(y)
    if np.any(nans):
        x = np.array(range(len(y)))
        y_nonan = np.interp(x, x[~nans], y[~nans])
    else:
        y_nonan = y

    y_filt = savgol_filter(y_nonan, window_size, order)

    # put any NaNs back
    if np.any(nans):
        y_filt[nans] = np.nan

    return y_filt


def run( soltab, axesToSmooth, size=[], mode='runningmedian', degree=1, replace=False, log=False, refAnt=''):
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

    refAnt : str, optional
        Reference antenna for phases. By default None.
    """

    import numpy as np
    from scipy.ndimage import generic_filter

    if refAnt == '': refAnt = None
    elif not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.warning('Reference antenna '+refAnt+' not found. Using: '+soltab.getAxisValues('ant')[1])
        refAnt = soltab.getAxisValues('ant')[1]

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
        vals = soltab.getValues(retAxesVals=False, refAnt=refAnt)
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
        for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToSmooth, weight=True, refAnt=refAnt):

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
