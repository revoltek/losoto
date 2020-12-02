#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading FLAG module.')

def _run_parser(soltab, parser, step):
    axesToFlag = parser.getarraystr( step, 'axesToFlag') # no default
    order = parser.getarrayint( step, 'order') # no default
    maxCycles = parser.getint( step, 'maxCycles', 5)
    maxRms = parser.getfloat( step, 'maxRms', 5.)
    maxRmsNoise = parser.getfloat( step, 'maxRmsNoise', 0.)
    fixRms = parser.getfloat( step, 'fixRms', 0.)
    fixRmsNoise = parser.getfloat( step, 'fixRmsNoise', 0.)
    windowNoise = parser.getint( step, 'windowNoise', 11)
    replace = parser.getbool( step, 'replace', False)
    preflagzeros = parser.getbool( step, 'preflagzeros', False)
    mode = parser.getstr( step, 'mode', 'smooth')
    refAnt = parser.getstr( step, 'refAnt', '')
    ncpu = parser.getint( '_global', 'ncpu', 0)

    parser.checkSpelling( step, soltab, ['axesToFlag', 'order', 'maxCycles', 'maxRms', 'maxRmsNoise', 'fixRms', 'fixRmsNoise', 'windowNoise', 'replace', 'preflagzeros', 'mode', 'refAnt'])
    return run( soltab, axesToFlag, order, maxCycles, maxRms, maxRmsNoise, fixRms, fixRmsNoise, windowNoise, replace, preflagzeros, mode, refAnt, ncpu )


def _flag(vals, weights, coord, solType, order, mode, preflagzeros, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRms, fixRmsNoise, replace, axesToFlag, selection, outQueue):

    import numpy as np
    import itertools
    from scipy.ndimage import generic_filter
    import scipy.interpolate

    def rolling_rms(a, window):
        """
        Return the rms for each element of the array calculated using the element inside a window,
        edges are mirrored
        """
        assert window % 2 == 1 # window must be odd
        a = np.pad(a, window/2, mode='reflect')
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        return np.sqrt(np.var(np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides), -1))

    def polyfit(x=None, y=None, z=None, w=None, order=None):
        """Two-dimensional polynomial fit.

        References:
            http://stackoverflow.com/questions/7997152/
            python-3d-polynomial-surface-fit-order-dependent/7997925#7997925

        x: array of coordinates
        y: array of coordinates (2d only)
        z: values of these coordinates (1d array as x and y)
        w: weights
        order: int for 1d, tuple of 2 values (x-order, y-order) for 2d

        Return: matrix of coefficents, for the 2d is ok to be fed into polynomial.polyval2d()
        """
        assert z.shape == w.shape
        if y is None:
            # note that np.polynomial.polynomial.polyfit want sigmas not variances
            m = np.polyfit(x, z, order, w=np.sqrt(w))

        else:
            from numpy.polynomial import polynomial
            x, y = np.meshgrid(np.asarray(x),np.asarray(y), indexing='ij')
            # TODO: unset x,y,z value if flagged
            z = z[(w != 0)]
            x = x[(w != 0)]
            y = y[(w != 0)]
            vander = polynomial.polyvander2d(x, y, order)
            vander = vander.reshape((-1,vander.shape[-1]))
            z = z.reshape((vander.shape[0],))
            m = np.linalg.lstsq(vander, z)[0]
            order = np.asarray(order)
            m = m.reshape(order+1) # matrix of coefficents

        return m

    def polyval(x=None, y=None, m=None):
        """Values to two-dimensional polynomial fit. Based uppon code
            provided by Joe Kington.
        """
        if y is None:
            poly = np.poly1d(m)
            return poly(x)

        else:
            x, y = np.meshgrid(x,y, indexing='ij')
            return np.polynomial.polynomial.polyval2d(x, y, m)


    def outlier_rej(vals, weights, axes, order=5, mode='smooth', max_ncycles=3, max_rms=3., max_rms_noise=0., window_noise=11., fix_rms=0., fix_rms_noise=0., replace=False):
        """
        Reject outliers using a running median
        val = the array (avg must be 0)
        weights = the weights to convert into flags
        axes = array with axes values (1d or 2d)
        order = "see polyfit()"
        max_ncycles = maximum number of cycles
        max_rms = number of rms times for outlier flagging
        max_rms_noise = cut on the rms of the rmss
        window_noise = window used to calculate the rmss to detect noise
        replace = instead of flag it, replace the data point with the smoothed one

        return: flags array and final rms
        """

        # renormalize axes to have decent numbers
        if len(axes) == 1:
            axes[0] -= axes[0][0]
            axes[0] /= axes[0][1] - axes[0][0]
        elif len(axes) == 2:
            axes[0] -= axes[0][0]
            axes[0] /= (axes[0][1]-axes[0][0])
            axes[1] -= axes[1][0]
            axes[1] /= (axes[1][1]-axes[1][0])
        else:
            logging.error('FLAG operation can flag only along 1 or 2 axes. Given axes: '+str(axesToFlag))
            return

        # test artificial data
        #axes[0] = np.array(range(len(axes[0])))*24414.0625
        #axes[1] = np.array(range(len(axes[1])))*5
        #vals = np.empty( shape=[len(axes[0]), len(axes[1])] )
        #for i, xi in enumerate(axes[0]):
        #    for j, yj in enumerate(axes[1]):
        #        vals[i,j] = 10*xi + yj**2
        #weights = np.ones_like(vals)

        if replace:
            orig_weights = np.copy(weights)

        for i in range(max_ncycles):

            # all is flagged? break
            if (weights == 0).all():
                rms = 0.
                break

            if mode == 'smooth':
                vals_smooth = np.copy(vals)
                np.putmask(vals_smooth, weights==0, np.nan)
                # speedup: if all data are used then just do a median and don't call the filter
                if all(o == 0 for o in order):
                    vals_smooth = np.ones(vals_smooth.shape)*np.nanmedian(vals_smooth)
                else:
                    for i, o in enumerate(order): 
                        if o == 0: order[i] = vals_smooth.shape[i]
                    vals_smooth = generic_filter(vals_smooth, np.nanmedian, size=order, mode='constant', cval=np.nan)
                vals_detrend = vals - vals_smooth
            # TODO: should be rolling
            elif mode == 'poly':
                # get polynomia and values
                if len(axes) == 1:
                    fit_sol = polyfit(axes[0], z=vals, w=weights, order=order)
                    vals_detrend = vals - polyval(axes[0], m=fit_sol)
                elif len(axes) == 2:
                    fit_sol = polyfit(axes[0], axes[1], z=vals, w=weights, order=order)
                    vals_smooth = polyval(axes[0], axes[1], m=fit_sol)
                    vals_detrend = vals - vals_smooth
            # TODO: should be rolling
            elif mode == 'spline':
                # get spline
                if len(axes) == 1:
                    spline = scipy.interpolate.UnivariateSpline(axes[0], y=vals, w=weights, k=order[0])
                    vals_detrend = vals - spline(axes[0])
                elif len(axes) == 2:
                    x, y = np.meshgrid(axes[0], axes[1], indexing='ij')
                    # spline doesn't like w=0
                    z = vals[(weights != 0)].flatten()
                    x = x[(weights != 0)].flatten()
                    y = y[(weights != 0)].flatten()
                    w = weights[(weights != 0)].flatten()
                    spline = scipy.interpolate.SmoothBivariateSpline(x, y, z, w, kx=order[0], ky=order[1])
                    vals_smooth = spline(axes[0], axes[1])
                    vals_detrend = vals - vals_smooth

            # remove outliers
            if max_rms > 0 or fix_rms > 0:
                # median calc https://en.wikipedia.org/wiki/Median_absolute_deviation
                rms =  1.4826 * np.nanmedian( np.abs(vals_detrend[(weights != 0)]) )
                if np.isnan(rms): weights[:] = 0
                elif fix_rms > 0:
                    flags = abs(vals_detrend) > fix_rms
                    weights[ flags ] = 0
                else:
                    flags = abs(vals_detrend) > max_rms * rms
                    weights[ flags ] = 0

            # remove noisy regions of data
            if max_rms_noise > 0 or fix_rms_noise > 0:
                rmses = rolling_rms(vals_detrend, window_noise)
                rms =  1.4826 * np.nanmedian( abs(rmses) )

                # rejection
                if fix_rms_noise > 0:
                    flags = rmses > fix_rms_noise
                else:
                    flags = rmses > (max_rms_noise * rms)
                weights[ flags ] = 0

            # all is flagged? break
            if (weights == 0).all():
                rms == 0.
                break

            # no flags? break
            if (flags == False).all():
                break

        # replace (outlier) flagged values with smoothed ones
        if replace:
            logging.debug('Replacing %.2f%% of the data.' % np.sum(orig_weights != weights))
            vals[np.where(orig_weights != weights)] = vals_smooth[np.where(orig_weights != weights)]
            weights = orig_weights

        # plot 1d
        plot = False
        if plot:
            import matplotlib as mpl
            mpl.use("Agg")
            import matplotlib.pyplot as plt
            plt.plot(axes[1][weights[0] == 0], vals[0][weights[0] == 0], 'ro')
            plt.plot(axes[1], vals[0], 'k.')
            plt.plot(axes[1], vals_smooth[0], 'r.')
            plt.plot(axes[1], vals_detrend[0], 'g.')
            plt.savefig('test.png')
            sys.exit(1)

        return weights, vals, rms


    def percentFlagged(w):
        return 100.*(weights.size-np.count_nonzero(weights))/float(weights.size)
    ########################################

    # check if everything flagged
    if (weights == 0).all() == True:
        logging.debug('Percentage of data flagged/replaced (%s): already completely flagged' % (removeKeys(coord, axesToFlag)))
        outQueue.put([vals, weights, selection])
        return

    if preflagzeros:
        if solType == 'amplitude': np.putmask(weights, vals == 1, 0)
        else: np.putmask(weights, vals == 0, 0)

    flagCoord = []
    for axisToFlag in axesToFlag:
        flagCoord.append(coord[axisToFlag])

    initPercentFlag = percentFlagged(weights)

    # works in phase-space (assume no wraps), remove just the mean to prevent problems if the phase is constantly around +/-pi
    if solType == 'phase' or solType == 'scalarphase' or solType == 'rotation':
        # remove mean of vals
        mean = np.angle( np.sum( weights.flatten() * np.exp(1j*vals.flatten()) ) / ( vals.flatten().size * sum(weights.flatten()) ) )
        logging.debug('Working in phase-space, remove angular mean '+str(mean)+'.')
        vals = normalize_phase(vals - mean)
        weights, vals, rms = outlier_rej(vals, weights, flagCoord, order, mode, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRms, fixRmsNoise, replace)
        vals = normalize_phase(vals + mean)

    elif solType == 'amplitude':
        vals_good = (vals>0)
        vals[vals_good] = np.log10(vals[vals_good])
        weights, vals, rms = outlier_rej(vals, weights, flagCoord, order, mode, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRmsNoise, replace)
        vals[vals_good] = 10**vals[vals_good]

    else:
        weights, vals, rms = outlier_rej(vals, weights, flagCoord, order, mode, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRms, fixRmsNoise, replace)

    if percentFlagged(weights) == initPercentFlag:
        logging.debug('Percentage of data flagged/replaced (%s): %.3f -> None' % ((removeKeys(coord, axesToFlag), initPercentFlag)))
    else:
        logging.debug('Percentage of data flagged/replaced (%s): %.3f -> %.3f %% (rms: %.5f)' \
            % ((removeKeys(coord, axesToFlag), initPercentFlag, percentFlagged(weights), rms)))

    outQueue.put([vals, weights, selection])
#    return vals, weights, selection


def run( soltab, axesToFlag, order, maxCycles=5, maxRms=5., maxRmsNoise=0., fixRms=0., fixRmsNoise=0., windowNoise=11, replace=False, preflagzeros=False, mode='smooth', refAnt='', ncpu=0 ):
    """
    This operation for LoSoTo implement a flagging procedure
    WEIGHT: compliant

    Parameters
    ----------
    axesToFlag : array of str
        Axes along which to smooth+find outlier (e.g. ['time', 'freq']), max 2 values.

    order : array of int
        Order of the function fitted during detrending. Array must have same size of axesToFlag. If mode=smooth these are the window of the running median (0=all axis).

    maxCycles : int, optional
        Max number of independent flagging cycles, by default 5.

    maxRms : float, optional
        Rms to clip outliers, by default 5.

    maxRmsNoise : float, optional
        Do a running rms and then flag those regions that have a rms higher than MaxRmsNoise*rms_of_rmses, by default 0 (ignored).

    fixRms : float, optional
        Instead of calculating rms use this value, by default 0 (ignored).

    fixRmsNoise : float, optional
        Instead of calculating rms of the rmses use this value (it will not be multiplied by the MaxRmsNoise), by default 0 (ignored).

    windowNoise : int, optional
        Window size for the running rms, by default 11.

    replace : bool, optional
        Replace bad values with the interpolated ones, instead of flagging them. By default False.

    preflagzeros : bool, optional
        Flag zeros/ones (bad solutions in BBS/DPPP). They should be flagged at import time. By default False.

    mode: str, optional
        Detrending/fitting algorithm: smooth / poly / spline. By default smooth.

    refAnt : str, optional
        Reference antenna, by default None.

    ncpu : int, optional
        Number of cpu to use, by default all available.
    """

    logging.info("Flag on soltab: "+soltab.name)

    # input check
    if refAnt == '':
        refAnt = None

    mode = mode.lower()
    if mode not in ['smooth', 'poly', 'spline']:
        logging.error('Mode must be one of smooth, poly, or spline')
        return 1

    if axesToFlag == []:
        logging.error("Please specify axis to flag. It must be a single one.")
        return 1

    for axisToFlag in axesToFlag:
        if axisToFlag not in soltab.getAxesNames():
            logging.error('Axis \"'+axisToFlag+'\" not found.')
            return 1

    if len(axesToFlag) != len(order) and (len(axesToFlag) != 1 or len(axesToFlag) != 2):
        logging.error("AxesToFlag and order must be both 1 or 2 values.")
        return 1

    if len(order) >= 2: order = list(order)

    # start processes for multi-thread
    mpm = multiprocManager(ncpu, _flag)

    # reorder axesToFlag as axes in the table
    axesToFlag_orig = axesToFlag
    axesToFlag = [coord for coord in soltab.getAxesNames() if coord in axesToFlag]
    if axesToFlag_orig != axesToFlag: order = order[::-1] # reverse order if we changed axesToFlag

    solType = soltab.getType()

    # fill the queue (note that sf and sw cannot be put into a queue since they have file references)
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToFlag, weight=True, refAnt=refAnt):
        mpm.put([vals, weights, coord, solType, order, mode, preflagzeros, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRms, fixRmsNoise, replace, axesToFlag, selection])

    mpm.wait()

    for v, w, sel in mpm.get():
        if replace:
            # rewrite solutions (flagged values are overwritten)
            soltab.setValues(v, sel, weight=False)
        else:
            soltab.setValues(w, sel, weight=True)

    soltab.flush()
    soltab.addHistory('FLAG (over %s with %s sigma cut)' % (axesToFlag, maxRms))

    return 0
