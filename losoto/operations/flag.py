#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a flagging procedure
# WEIGHT: flag-only compliant

import logging
from losoto.operations_lib import *
import numpy as np
import itertools
from scipy.ndimage import generic_filter

logging.debug('Loading FLAG module.')

def flag(vals, weights, coord, solType, order, smooth, preflagzeros, maxCycles, maxRms, replace, axesToFlag, selection):
#def flag(vals, weights, coord, solType, order, smooth, preflagzeros, maxCycles, maxRms, replace, axesToFlag, selection, outQueue):

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
            print z.shape, x.shape, y.shape
            m = np.linalg.lstsq(vander, z)[0]
            order = np.asarray(order)
            m = m.reshape(order+1) # matrix of coefficents

        return m

    def  polyval(x=None, y=None, m=None):
        """Values to two-dimensional polynomial fit. Based uppon code 
            provided by Joe Kington.
        """
        if y is None:
            poly = np.poly1d(m)
            return poly(x)

        else:
            x, y = np.meshgrid(x,y, indexing='ij')
            return np.polynomial.polynomial.polyval2d(x, y, m)
 

#    def clean_noisy(data, times, window, max_rms):
#        """
#        calculate a running RMS and remove noisy data
#
#        window = in timestamps, sliding window dimension
#        max_rms = flag points in a region with rms larger than max_rms times the rms of rmses
#        
#        return: an array of data dimensions with flags
#        """
#        if len(data) == 0: return []
#        # loop over solution times
#        rmses = np.zeros(shape=data.shape, dtype=np.float)
#        for i, time in enumerate(times):
#    
#            # get data to smooth (values inside the time window)
#            data_array = data[ np.where( abs(times - time) <= window / 2. ) ]
#            rmses[i] = np.std(data_array)
#
#        rms =  1.4826 * np.median( abs(rmses) )
#        flags = abs(rmses) > max_rms * rms
# 
#        return flags
    
    def outlier_rej(vals, weights, axes, order=5, smooth=False, max_ncycles = 3, max_rms = 3., replace = False):
        """
        Reject outliers using a running median
        val = the array (avg must be 0)
        weights = the weights to convert into flags
        axes = array with axes values (1d or 2d)
        order = "see polyfit()"
        max_ncycles = maximum number of cycles
        max_rms = number of rms times for outlier flagging
        replace = instead of flag it, replace the data point with the smoothed one
    
        return: flags array and final rms
        """
      
        # renormalize axes to have decent numbers
        if len(axes) == 1:
            axes[0] -= axes[0][0]
            axes[0] /= axes[0][1] - axes[0][0]
            order = order[0]
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
                rms == 0.
                break

            if smooth:
                    vals_smooth = np.copy(vals)
                    np.putmask(vals_smooth, weights==0, np.nan)
                    vals_smooth = generic_filter(vals_smooth, np.nanmedian, size=order)
                    vals_detrend = vals - vals_smooth
            else:
                # get polynomia and values
                if len(axes) == 1: 
                    fit_sol = polyfit(axes[0], z=vals, w=weights, order=order)
                    vals_detrend = vals - polyval(axes[0], m=fit_sol)
                elif len(axes) == 2: 
                    fit_sol = polyfit(axes[0], axes[1], z=vals, w=weights, order=order)
                    vals_detrend = vals - polyval(axes[0], axes[1], m=fit_sol)

            # remove noisy regions of data
            #flag_noisy = clean_noisy(vals_detrend, time[ s ], window, max_rms_noise)
            #print 'noise flagging', float(sum(flag_noisy))/len(flag_noisy)
            #flags[ s ] = flag_noisy # add flags g (shape s=True) to global flags
            #s[ s ] = ~flag_noisy # new refined selection
            #vals_detrend = vals[ s ] - vals_smoothed[ ~flag_noisy ] # keep only vals satisfying s and g

            # median calc
            rms =  1.4826 * np.median( np.abs(vals_detrend[(weights != 0)]) )
    
            # rejection  
            flags = abs(vals_detrend) > max_rms * rms
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
        #import matplotlib.pyplot as plt
        #plt.plot(axes[0], vals, 'k.')
        #plt.plot(axes[0][weights == 0], vals[weights == 0], 'ro')
        #plt.plot(axes[0], vals_smooth, 'r-')
        #plt.savefig('test.png')
        #sys.exit(1)

        # plot 2d
        #import matplotlib.pyplot as plt
        #plt.imshow(vals.T, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=1./5)
        #plt.colorbar()
        #plt.savefig('test2d.png')
        #plt.clf()
        #plt.imshow(polyval(axes[0], axes[1], m=fit_sol).T, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=1./5)
        #plt.imshow(vals_smooth.T, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=1./5)
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
        return vals, weights, selection
#        outQueue.put([vals, weights, selection])
#        return

    if preflagzeros:
        if solType == 'amplitude': weights[np.where(vals == 1)] = 0
        else: np.putmask(weights, vals == 0, 0)

    flagCoord = []
    for axisToFlag in axesToFlag:
        flagCoord.append(coord[axisToFlag])

    initPercentFlag = percentFlagged(weights)

    # if phase, then convert to real/imag, run the flagger on those, and convert back to pahses
    # best way to avoid unwrapping
    if solType == 'phase' or solType == 'scalarphase' or solType == 'rotation':
        re = 1. * np.cos(vals)
        im = 1. * np.sin(vals)
        weights_re, re, rms_re = outlier_rej(re, weights, flagCoord, order, smooth, maxCycles, maxRms, replace)
        weights_im, im, rms_im = outlier_rej(im, weights, flagCoord, order, smooth, maxCycles, maxRms, replace)
        vals = np.arctan2(im, re)
        np.putmask(weights, weights_re == 0, 0)
        np.putmask(weights, weights_im == 0, 0)
        rms = np.sqrt(rms_re**2 + rms_im**2)

    elif solType == 'amplitude':
        weights, vals, rms = outlier_rej(np.log10(vals), weights, flagCoord, order, smooth, maxCycles, maxRms, replace)
        vals == 10**vals

    else:
        weights, vals, rms = outlier_rej(vals, weights, flagCoord, order, smooth, maxCycles, maxRms, replace)
    
    if percentFlagged(weights) == initPercentFlag:
        logging.debug('Percentage of data flagged/replaced (%s): %.3f -> None' % (removeKeys(coord, axisToFlag), initPercentFlag))
    else: 
        logging.debug('Percentage of data flagged/replaced (%s): %.3f -> %.3f %% (rms: %.5f)' \
            % (removeKeys(coord, axisToFlag), initPercentFlag, percentFlagged(weights), rms))

#    outQueue.put([vals, flags, selection])
    return vals, weights, selection
        
            
def run( step, parset, H ):

    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    #check_parset('Axes','MaxCycles','MaxRms','Order','Replace','PreFlagZeros')
    axesToFlag = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), 'time' )
    maxCycles = parset.getInt('.'.join(["LoSoTo.Steps", step, "MaxCycles"]), 5 )
    maxRms = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRms"]), 5. )
    order = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "Order"]), 3 )
    replace = parset.getBool('.'.join(["LoSoTo.Steps", step, "Replace"]), False )
    preflagzeros = parset.getBool('.'.join(["LoSoTo.Steps", step, "PreFlagZeros"]), False )
    smooth = parset.getBool('.'.join(["LoSoTo.Steps", step, "Smooth"]), False )
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 1 )

    if axesToFlag == []:
        logging.error("Please specify axis to flag. It must be a single one.")
        return 1

    if len(axesToFlag) != len(order):
        logging.error("AxesToFlag and order must be both 1 or 2 values.")
        return 1

    if len(order) == 1: order = order[0]
    elif len(order) == 2: order = tuple(order)


    # start processes for multi-thread
#    mpm = multiprocManager(ncpu, flag)

    for soltab in openSoltabs( H, soltabs ):

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
                return 1

        # reorder axesToFlag as axes in the table
        axesToFlag = [axisToFlag for (coord, axisToFlag) in zip(sf.getAxesNames(),axesToFlag)]
        if type(order) is int: order = [order]
        order = [order for (coord, order) in zip(sf.getAxesNames(),order)]

        solType = sf.getType()

        # fill the queue (note that sf and sw cannot be put into a queue since they have file references)
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=axesToFlag, weight=True):
            #mpm.put([vals, weights, coord, solType, order, smooth, preflagzeros, maxCycles, maxRms, replace, axesToFlag, selection])
            v, w, sel = flag(vals, weights, coord, solType, order, smooth, preflagzeros, maxCycles, maxRms, replace, axesToFlag, selection)

#        mpm.wait()
        
#        for v, w, sel in mpm.get():
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
