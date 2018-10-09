#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

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
    maxFlaggedFraction = parser.getfloat( step, 'maxFlaggedFraction', 1.0)
    maxStddev = parser.getfloat( step, 'maxStddev', 0.01)
    replace = parser.getbool( step, 'replace', False)
    preflagzeros = parser.getbool( step, 'preflagzeros', False)
    mode = parser.getstr( step, 'mode', 'smooth')
    telescope = parser.getstr( step, 'telescope', 'lofar')
    refAnt = parser.getstr( step, 'refAnt', '')
    ncpu = parser.getint( '_global', 'ncpu', 0)

    parser.checkSpelling( step, soltab, ['axesToFlag', 'order', 'maxCycles', 'maxRms', 'maxRmsNoise', 'fixRms', 'fixRmsNoise', 'windowNoise', 'maxFlaggedFraction', 'maxStddev', 'replace', 'preflagzeros', 'mode', 'telescope', 'refAnt'])
    return run( soltab, axesToFlag, order, maxCycles, maxRms, maxRmsNoise, fixRms, fixRmsNoise, windowNoise, maxFlaggedFraction, maxStddev, replace, preflagzeros, mode, telescope, refAnt, ncpu )


def _flag(vals, weights, coord, solType, order, mode, preflagzeros, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRms, fixRmsNoise, replace, axesToFlag, maxFlaggedFraction, selection, outQueue):

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

        for i in xrange(max_ncycles):

            # all is flagged? break
            if (weights == 0).all():
                rms = 0.
                break

            if mode == 'smooth':
                vals_smooth = np.copy(vals)
                np.putmask(vals_smooth, weights==0, np.nan)
                if all(o == 0 for o in order): order = vals_smooth.shape
                vals_smooth = generic_filter(vals_smooth, np.nanmedian, size=order)
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
            plt.plot(axes[0][weights == 0], vals[weights == 0], 'ro')
            plt.plot(axes[0], vals, 'k.')
            plt.plot(axes[0], vals_smooth, 'r.')
            #plt.plot(axes[0], vals_detrend, 'g.')
            plt.savefig('test.png')
            #sys.exit(1)

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
    elif percentFlagged(weights)/100.0 >= maxFlaggedFraction:
        logging.debug('Percentage of data flagged (%s): %.3f -> %.3f %% (rms: %.5f)' \
            % ((removeKeys(coord, axesToFlag), initPercentFlag, percentFlagged(weights), rms)))
        logging.debug('Flagged percentage exceeds maxFlaggedFraction. Flagging all')
        weights[:] = 0
    else:
        logging.debug('Percentage of data flagged/replaced (%s): %.3f -> %.3f %% (rms: %.5f)' \
            % ((removeKeys(coord, axesToFlag), initPercentFlag, percentFlagged(weights), rms)))

    outQueue.put([vals, weights, selection])
#    return vals, weights, selection


def _flag_bandpass(freqs, amps, weights, telescope, nSigma, maxFlaggedFraction, maxStddev,
                     plot, s, outQueue):
    """
    Flags bad amplitude solutions relative to median bandpass (in log space) by setting
    the corresponding weights to 0.0

    Note: A median over the time axis is done before flagging, so the flags are not time-
    dependent

    Parameters
    ----------
    freqs : array
        Array of frequencies

    amps : array
        Array of amplitudes as [time, ant, freq, pol]

    weights : array
        Array of weights as [time, ant, freq, pol]

    telescope : str, optional
        Specifies the telescope for the bandpass model

    nSigma : float
        Number of sigma for flagging. Amplitudes outside of nSigma*stddev are flagged

    maxFlaggedFraction : float
        Maximum allowable fraction of flagged frequencies. Stations with higher fractions
        will be completely flagged

    maxStddev : float
        Maximum allowable standard deviation

    plot : bool
        If True, the bandpass with flags and best-fit line is plotted for each station

    s : int
        Station index

    Returns
    -------
    indx, weights : int, array
        Station index, modified weights array
    """
    def _B(x, k, i, t, extrap, invert):
        if k == 0:
            if extrap:
                if invert:
                    return -1.0
                else:
                    return 1.0
            else:
                return 1.0 if t[i] <= x < t[i+1] else 0.0
        if t[i+k] == t[i]:
           c1 = 0.0
        else:
           c1 = (x - t[i])/(t[i+k] - t[i]) * _B(x, k-1, i, t, extrap, invert)
        if t[i+k+1] == t[i+1]:
           c2 = 0.0
        else:
           c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * _B(x, k-1, i+1, t, extrap, invert)
        return c1 + c2


    def _bspline(x, t, c, k):
        n = len(t) - k - 1
        assert (n >= k+1) and (len(c) >= n)
        invert = False
        extrap = [False] * n
        if x >= t[n]:
            extrap[-1] = True
        elif x < t[k]:
            extrap[0] = True
            invert = False
        return sum(c[i] * _B(x, k, i, t, e, invert) for i, e in zip(range(n), extrap))


    def _bandpass_LBA(freq, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13):
        """
        Defines the functional form of the LBA bandpass in terms of splines of degree 3

        The spline fit was done using LSQUnivariateSpline() on the median bandpass between
        30 MHz and 78 MHz. The knots were set by hand to acheive a good fit with a
        minimum number of parameters.

        Parameters
        ----------
        freq : array
            Array of frequencies

        c1-c13 : float
            Spline coefficients

        Returns
        -------
        bandpass : list
            List of bandpass values as function of frequency
        """
        knots = np.array([30003357.0, 30003357.0, 30003357.0, 30003357.0, 40000000.0,
                          50000000.0, 55000000.0, 56000000.0, 60000000.0, 62000000.0,
                          63000000.0, 64000000.0, 70000000.0, 77610779.0, 77610779.0,
                          77610779.0, 77610779.0])
        coeffs = np.array([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13])
        return [_bspline(f, knots, coeffs, 3) for f in freq]


    def _bandpass_HBA_low(freq, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10):
        """
        Defines the functional form of the HBA-low bandpass in terms of splines of degree
        3

        The spline fit was done using LSQUnivariateSpline() on the median bandpass between
        120 MHz and 188 MHz. The knots were set by hand to acheive a good fit with a
        minimum number of parameters.

        Parameters
        ----------
        freq : array
            Array of frequencies

        c1-c10 : float
            Spline coefficients

        Returns
        -------
        bandpass : list
            List of bandpass values as function of frequency
        """
        knots = np.array([1.15e+08, 1.15e+08, 1.15e+08, 1.15e+08,
                          1.30e+08, 1.38e+08, 1.48e+08, 1.60e+08,
                          1.68e+08, 1.78e+08, 1.90e+08, 1.90e+08,
                          1.9e+08, 1.9e+08])
        coeffs = np.array([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10])
        return [_bspline(f, knots, coeffs, 3) for f in freq]


    def _fit_bandpass(freq, logamp, sigma, band, do_fit=True):
        """
        Fits amplitudes with one of the bandpass functions

        The initial coefficients were determined from a LSQUnivariateSpline() fit on the
        median bandpass of the appropriate band. The allowable fitting ranges were set by
        hand through testing on a number of observations (to allow the bandpass function
        to adjust for the differences between stations but not to fit to RFI, etc.).

        Parameters
        ----------
        freq : array
            Array of frequencies

        amps : array
            Array of log10(amplitudes)

        sigma : array
            Array of sigma (1/weights**2)

        band : str
            Band name ('hba_low', etc.)

        do_fit : bool, optional
            If True, the fitting is done. If False, the unmodified model bandpass is
            returned

        Returns
        -------
        fit_parms, bandpass : list, list
            List of best-fit parameters, List of bandpass values as function of frequency
        """
        from scipy.optimize import curve_fit

        if band.lower() == 'hba_low':
            bandpass_function = _bandpass_HBA_low
            init_coeffs = np.array([-0.01460369, 0.05062699, 0.02827004, 0.03738518,
                                    -0.05729109, 0.02303295, -0.03550487, -0.0803113,
                                    -0.2394929, -0.358301])
            bounds_deltas_lower = [0.06, 0.05, 0.04, 0.04, 0.04, 0.04, 0.1, 0.1, 0.2, 0.5]
            bounds_deltas_upper = [0.06, 0.1, 0.1, 0.1, 0.04, 0.04, 0.04, 0.04, 0.05, 0.06]
        elif band.lower() == 'lba':
            bandpass_function = _bandpass_LBA
            init_coeffs = np.array([-0.22654016, -0.1950495, -0.07763014, 0.10002095,
                                    0.32797671, 0.46900048, 0.47155583, 0.31945897,
                                    0.29072278, 0.08064795, -0.15761538, -0.36020451,
                                    -0.51163338])
            bounds_deltas_lower = [0.25, 0.2, 0.05, 0.05, 0.05, 0.05, 0.08, 0.05, 0.08, 0.15,
                                   0.15, 0.15, 0.15]
            bounds_deltas_upper = [0.25, 0.2, 0.05, 0.05, 0.05, 0.05, 0.08, 0.05, 0.08, 0.15,
                                   0.15, 0.15, 0.15]
        else:
            print('The "{}" band is not supported'.format(band))
            sys.exit(1)

        if do_fit:
            lower = [c - b for c, b in zip(init_coeffs, bounds_deltas_lower)]
            upper = [c + b for c, b in zip(init_coeffs, bounds_deltas_upper)]
            param_bounds = (lower, upper)
            try:
                popt, pcov = curve_fit(bandpass_function, freq, logamp, sigma=sigma,
                                       bounds=param_bounds, method='dogbox')
                return popt, bandpass_function(freq, *tuple(popt))
            except RuntimeError:
                logging.error('Fitting failed.' )
                return None, bandpass_function(freq, *tuple(init_coeffs))
        else:
            return None, bandpass_function(freq, *tuple(init_coeffs))

    # Check that telescope is supported
    if telescope.lower() == 'lofar':
        # Determine which band we're in
        if np.median(freqs) < 180e6 and np.median(freqs) > 110e6:
            band = 'hba_low'
            median_min = 50.0
            median_max = 200.0
        elif np.median(freqs) < 90e6:
            band = 'lba'
            median_min = 50.0
            median_max = 200.0
        else:
            print('The median frequency of {} Hz is outside of any supported LOFAR band '
                  '(LBA and HBA-low)'.format(np.median(freqs)))
            sys.exit(1)
    else:
       logging.error("Only telescope = 'lofar' is currently supported for bandpass mode.")
       outQueue.put([s, weights])
       return 1

    # Skip fully flagged stations
    if np.all(weights == 0.0):
        outQueue.put([s, weights])
        return

    # Build arrays for fitting
    flagged = np.where(np.logical_or(weights == 0.0, amps == 0.0))
    amps_flagged = amps.copy()
    amps_flagged[flagged] = np.nan
    sigma = weights.copy()
    sigma[flagged] = 1.0
    sigma = np.sqrt(1.0 / sigma)
    sigma[flagged] = 1e8

    # Iterate over polarizations
    npols = amps.shape[2] # number of polarizations
    for pol in range(npols):
        # take median over time and divide out the median offset
        with np.warnings.catch_warnings():
            # filter NaN warnings -- we deal with NaNs below
            np.warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
            amps_div = np.nanmedian(amps_flagged[:, :, pol], axis=0)
            median_val = np.nanmedian(amps_div)
        amps_div /= median_val
        sigma_div = np.median(sigma[:, :, pol], axis=0)
        median_flagged = np.where(np.isnan(amps_div))
        amps_div[median_flagged] = 1.0
        sigma_div[median_flagged] = 1e8
        median_flagged = np.where(amps_div <= 0.0)
        amps_div[median_flagged] = 1.0
        sigma_div[median_flagged] = 1e8
        sigma_orig = sigma_div.copy()

        # Before doing the fitting, flag any solutions that deviate from the model bandpass by
        # a large factor to avoid biasing the first fit
        _, bp_sp = _fit_bandpass(freqs, np.log10(amps_div), sigma_div, band, do_fit=False)
        bad = np.where(np.abs(bp_sp - np.log10(amps_div)) > 0.2)
        sigma_div[bad] = 1e8

        # Iteratively fit and flag
        maxiter = 5
        niter = 0
        nflag = 0
        nflag_prev = -1
        while nflag != nflag_prev and niter < maxiter:
            p, bp_sp = _fit_bandpass(freqs, np.log10(amps_div), sigma_div, band)
            stdev_all = np.sqrt(np.average((bp_sp-np.log10(amps_div))**2, weights=(1/sigma_div)**2))
            stdev = min(maxStddev, stdev_all)
            bad = np.where(np.abs(bp_sp - np.log10(amps_div)) > nSigma*stdev)
            nflag = len(bad[0])
            if nflag == 0:
                break
            if niter > 0:
                nflag_prev = nflag
            sigma_div = sigma_orig.copy()  # reset flags to original ones
            sigma_div[bad] = 1e8
            niter += 1

        if plot:
            import matplotlib.pyplot as plt
            plt.plot(freqs, bp_sp, 'g-', lw=3)
            plt.plot(freqs, np.log10(amps_div), 'o', c='g')
            plt.plot(freqs[bad], np.log10(amps_div)[bad], 'o', c='r')
            plt.show()

        # Check whether entire station is bad (high stdev or high flagged fraction). If
        # so, flag all frequencies and polarizations
        if stdev_all > maxStddev * 5.0:
            # Station has high stddev relative to median bandpass
            logging.info('Flagged station {0} (pol {1}) due to high stddev '
                  '({2})'.format(s, pol, stdev_all))
            weights[:, :, pol] = 0.0
        elif float(len(bad[0]))/float(len(freqs)) > maxFlaggedFraction:
            # Station has high fraction of flagged solutions
            logging.info('Flagged station {0} (pol {1}) due to high flagged fraction '
                  '({2})'.format(s, pol, float(len(bad[0]))/float(len(freqs))))
            weights[:, :, pol] = 0.0
        elif median_val < median_min or median_val > median_max:
            # Station has extreme median value
            logging.info('Flagged station {0} (pol {1}) due to extreme median value '
                  '({2})'.format(s, pol, median_val))
            weights[:, :, pol] = 0.0
        else:
            # Station is OK; flag solutions with high sigma values
            flagged = np.where(sigma_div > 1e3)
            nflagged_orig = len(np.where(weights[:, :, pol] == 0.0)[0])
            weights[:, flagged[0], pol] = 0.0
            nflagged_new = len(np.where(weights[:, :, pol] == 0.0)[0])
            prcnt = float(nflagged_new - nflagged_orig) / float(np.product(weights.shape[:-1])) * 100.0
            logging.info('Flagged {0}% of solutions for station {1} (pol {2})'.format(prcnt, s, pol))

    outQueue.put([s, weights])


def run( soltab, axesToFlag, order, maxCycles=5, maxRms=5., maxRmsNoise=0., fixRms=0., fixRmsNoise=0., windowNoise=11, maxFlaggedFraction=1.0, maxStddev=0.01, replace=False, preflagzeros=False, mode='smooth', telescope='lofar', refAnt='', ncpu=0 ):
    """
    This operation for LoSoTo implement a flagging procedure
    WEIGHT: compliant

    Parameters
    ----------
    axesToFlag : array of str
        Axes along which to smooth+find outlier (e.g. ['time', 'freq']), max 2 values. Must be ['freq'] if mode=bandpass.

    order : array of int
        Order of the function fitted during detrending. Array must have same size of axesToFlag. If mode=smooth these are the window of the running median (0=all axis). If mode=bandpass this value is ignored.

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

    maxFlaggedFraction : float, optional
        If mode=bandpass, this sets the maximum allowable fraction of flagged solutions
        above which the entire station is flagged.

    maxStddev : float, optional
        If mode=bandpass, this sets the maximum allowable standard deviation considered when outlier clipping is done

    replace : bool, optional
        Replace bad values with the interpolated ones, instead of flagging them. By default False.

    preflagzeros : bool, optional
        Flag zeros/ones (bad solutions in BBS/DPPP). They should be flagged at import time. By default False.

    mode: str, optional
        Detrending/fitting algorithm: smooth / poly / spline / bandpass. By default smooth.

    telescope : str, optional
        Specifies the telescope if mode = 'bandpass'.

    refAnt : str, optional
        Reference antenna, by default None.

    ncpu : int, optional
        Number of cpu to use, by default all available.
    """

    logging.info("Flag on soltab: "+soltab.name)

    # input check
    if ncpu == 0:
        import multiprocessing
        ncpu = multiprocessing.cpu_count()

    if refAnt == '':
        refAnt = None

    mode = mode.lower()
    if mode not in ['smooth', 'poly', 'spline', 'bandpass']:
        logging.error('Mode must be one of smooth, poly, spline, or bandpass')
        return 1

    if axesToFlag == []:
        logging.error("Please specify axis to flag. It must be a single one.")
        return 1

    for axisToFlag in axesToFlag:
        if axisToFlag not in soltab.getAxesNames():
            logging.error('Axis \"'+axisToFlag+'\" not found.')
            return 1

    if mode == 'bandpass':
        solType = soltab.getType()
        if solType != 'amplitude':
           logging.error("Soltab must be of type amplitude for bandpass mode.")
           return 1
        if replace:
           logging.warning("replace = True is not currently supported in bandpass mode. Ignoring")
        if axesToFlag != ['freq']:
            logging.error('For bandpass mode, the axis to flag must be freq.')
            return 1

        # Axis order must be [time, ant, freq, pol], so reorder if necessary
        axis_names = soltab.getAxesNames()
        if ('freq' not in axis_names or 'pol' not in axis_names or
            'time' not in axis_names or 'ant' not in axis_names):
           logging.error("Currently, bandpass mode requires the following axes: "
                         "freq, pol, time, and ant.")
           return 1
        freq_ind = axis_names.index('freq')
        pol_ind = axis_names.index('pol')
        time_ind = axis_names.index('time')
        ant_ind = axis_names.index('ant')
        amps_arraytmp = soltab.val[:].transpose([time_ind, ant_ind, freq_ind, pol_ind])
        weights_arraytmp = soltab.weight[:].transpose([time_ind, ant_ind, freq_ind, pol_ind])
        if preflagzeros:
            flagged = np.where(amplitude_arraytmp == 1.0)
            weights_arraytmp[flagged] = 0.0

        # fill the queue
        mpm = multiprocManager(ncpu, _flag_bandpass)
        for s in range(len(soltab.ant)):
            mpm.put([soltab.freq[:], amps_arraytmp[:, s, :, :], weights_arraytmp[:, s, :, :],
                     telescope, maxRms, maxFlaggedFraction, maxStddev, False, s])
        mpm.wait()

        # Write new weights
        for (s, w) in mpm.get():
            weights_arraytmp[:, s, :, :] = w
        weights_array = weights_arraytmp.transpose([time_ind, ant_ind, freq_ind, pol_ind])
        soltab.setValues(weights_array, weight=True)
        soltab.addHistory('FLAG (mode=bandpass, telescope={0}, maxRms={1}, '
                          'maxFlaggedFraction={2}, maxStddev={3}'.format(telescope, maxRms,
                          maxFlaggedFraction, maxStddev))
    else:
        if len(axesToFlag) != len(order) and (len(axesToFlag) != 1 or len(axesToFlag) != 2):
            logging.error("AxesToFlag and order must be both 1 or 2 values.")
            return 1

        if len(order) == 2: order = tuple(order)

        # start processes for multi-thread
        mpm = multiprocManager(ncpu, _flag)

        # reorder axesToFlag as axes in the table
        axesToFlag_orig = axesToFlag
        axesToFlag = [coord for coord in soltab.getAxesNames() if coord in axesToFlag]
        if axesToFlag_orig != axesToFlag: order = order[::-1] # reverse order if we changed axesToFlag

        solType = soltab.getType()

        # fill the queue (note that sf and sw cannot be put into a queue since they have file references)
        for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToFlag, weight=True, reference=refAnt):
            mpm.put([vals, weights, coord, solType, order, mode, preflagzeros, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRms, fixRmsNoise, replace, axesToFlag, maxFlaggedFraction, selection])

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
