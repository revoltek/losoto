#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a flagging procedure
# WEIGHT: flag-only compliant

import logging
from losoto.operations_lib import *
import numpy as np
import itertools
from scipy.ndimage import generic_filter
import scipy.interpolate

logging.debug('Loading FLAG module.')

def flag(vals, weights, coord, solType, order, mode, preflagzeros, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRmsNoise, replace, axesToFlag, selection, outQueue):
    
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

    def normalize(phase):
        """
        Normalize phase to the range [-pi, pi].
        """
        # Convert to range [-2*pi, 2*pi].
        out = np.fmod(phase, 2.0 * np.pi)
        # Remove nans
        np.putmask(out, out!=out, 0)
        # Convert to range [-pi, pi]
        out[out < -np.pi] += 2.0 * np.pi
        out[out > np.pi] -= 2.0 * np.pi
        return out


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
 

    def outlier_rej(vals, weights, axes, order=5, mode='smooth', max_ncycles=3, max_rms=3., max_rms_noise=0., window_noise=11., fix_rms_noise=0., replace=False):
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
                    vals_detrend = vals - polyval(axes[0], axes[1], m=fit_sol)
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
                    vals_detrend = vals - spline(axes[0], axes[1])

            # remove outliers
            if max_rms > 0:
                # median calc https://en.wikipedia.org/wiki/Median_absolute_deviation
                rms =  1.4826 * np.nanmedian( np.abs(vals_detrend[(weights != 0)]) )
                if np.isnan(rms): weights[:] = 0
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
            vals[orig_weights != weights] = vals_detrended
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

            # plot 2d
            #import matplotlib as mpl
            #mpl.use("Agg")
            #import matplotlib.pyplot as plt
            #plt.imshow(vals.T, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=1./5)
            #plt.colorbar()
            #plt.savefig('test2d.png')
            #plt.clf()
            ##plt.imshow(polyval(axes[0], axes[1], m=fit_sol).T, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=1./5)
            ##plt.imshow(vals_smooth.T, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=1./5)
            #plt.imshow(spline(axes[0],axes[1]).T, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=1./5)
            #plt.colorbar()
            #plt.savefig('test2d-smooth.png')
            #plt.clf()
            #plt.imshow(vals_detrend.T, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=1/5.)
            #plt.colorbar()
            #plt.savefig('test2d-detrend.png')
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
        vals = normalize(vals - mean)
        weights, vals, rms = outlier_rej(vals, weights, flagCoord, order, mode, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRmsNoise, replace)
        vals = normalize(vals + mean)

    elif solType == 'amplitude':
        weights, vals, rms = outlier_rej(np.log10(vals), weights, flagCoord, order, mode, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRmsNoise, replace)
        vals = 10**vals

    else:
        weights, vals, rms = outlier_rej(vals, weights, flagCoord, order, mode, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRmsNoise, replace)
    
    clean_coord = {key: coord[key] for key in coord if key not in axesToFlag}
    if percentFlagged(weights) == initPercentFlag:
        logging.debug('Percentage of data flagged/replaced (%s): %.3f -> None' % (clean_coord, initPercentFlag))
    else: 
        logging.debug('Percentage of data flagged/replaced (%s): %.3f -> %.3f %% (rms: %.5f)' \
            % (clean_coord, initPercentFlag, percentFlagged(weights), rms))

    outQueue.put([vals, weights, selection])
#    return vals, weights, selection
        
            
def run( step, parset, H ):

    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    #check_parset('Axes','MaxCycles','MaxRms','Order','Replace','PreFlagZeros')
    axesToFlag = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), 'time' )
    maxCycles = parset.getInt('.'.join(["LoSoTo.Steps", step, "MaxCycles"]), 5 )
    maxRms = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRms"]), 5. )
    maxRmsNoise = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRmsNoise"]), 0. )
    fixRmsNoise = parset.getFloat('.'.join(["LoSoTo.Steps", step, "FixRmsNoise"]), 0. )
    windowNoise = parset.getInt('.'.join(["LoSoTo.Steps", step, "WindowNoise"]), 11 )
    order = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "Order"]), 3 )
    replace = parset.getBool('.'.join(["LoSoTo.Steps", step, "Replace"]), False )
    preflagzeros = parset.getBool('.'.join(["LoSoTo.Steps", step, "PreFlagZeros"]), False )
    mode = parset.getString('.'.join(["LoSoTo.Steps", step, "Mode"]), 'smooth' )
    ref = parset.getString('.'.join(["LoSoTo.Steps", step, "Reference"]), '' )
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 0 )
    if ncpu == 0:
        import multiprocessing
        ncpu = multiprocessing.cpu_count()

    if ref == '': ref = None

    if axesToFlag == []:
        logging.error("Please specify axis to flag. It must be a single one.")
        return 1

    if len(axesToFlag) != len(order):
        logging.error("AxesToFlag and order must be both 1 or 2 values.")
        return 1

    if len(order) == 2: order = tuple(order)

    mode = mode.lower()
    if mode != 'smooth' and mode != 'poly' and mode != 'spline':
        logging.error('Mode must be smooth, poly or spline')
        return 1

    for soltab in openSoltabs( H, soltabs ):

        # start processes for multi-thread
        mpm = multiprocManager(ncpu, flag)

        logging.info("Flagging soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab, useCache=True) # remember to flush!

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        for axisToFlag in axesToFlag:
            if axisToFlag not in sf.getAxesNames():
                logging.error('Axis \"'+axis+'\" not found.')
                mpm.wait()
                return 1

        # reorder axesToFlag as axes in the table
        axesToFlag_orig = axesToFlag
        axesToFlag = [coord for coord in sf.getAxesNames() if coord in axesToFlag]
        if axesToFlag_orig != axesToFlag: order = order[::-1] # reverse order if we changed axesToFlag

        solType = sf.getType()

        # fill the queue (note that sf and sw cannot be put into a queue since they have file references)
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=axesToFlag, weight=True, reference=ref):
            mpm.put([vals, weights, coord, solType, order, mode, preflagzeros, maxCycles, maxRms, maxRmsNoise, windowNoise, fixRmsNoise, replace, axesToFlag, selection])

        mpm.wait()
        
        for v, w, sel in mpm.get():
            sw.selection = sel
            if replace:
                # rewrite solutions (flagged values are overwritten)
                sw.setValues(v, weight=False)
            else:
                sw.setValues(w, weight=True)
        
        sw.flush()
        sw.addHistory('FLAG (over %s with %s sigma cut)' % (axesToFlag, maxRms))

        del sw
        del sf
        del soltab

    return 0
