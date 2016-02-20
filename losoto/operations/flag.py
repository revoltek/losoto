#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a flagging procedure
# WEIGHT: flag-only compliant

import logging
from losoto.operations_lib import *
import numpy as np
import itertools

logging.debug('Loading FLAG module.')

def flag(vals, weights, coord, solType, order, preflagzeros, maxCycles, maxRms, maxRmsNoise, replace, axesToFlag, selection):
#def flag(vals, weights, coord, solType, order, preflagzeros, maxCycles, maxRms, maxRmsNoise, replace, axesToFlag, selection, outQueue):

    def polyfit(x=None, y=None, z=None, w=None, order=3, typ='1d'):
        """Two-dimensional polynomial fit. Based uppon code provided by 
        Joe Kington.
    
        References:
            http://stackoverflow.com/questions/7997152/
            python-3d-polynomial-surface-fit-order-dependent/7997925#7997925

        """
        assert z.shape == weights.shape
        if typ == '1d':
            # note that np.polynomial.polynomial.polyfit want sigmas not variances
            m = np.polyfit(x, z, order, w=np.sqrt(w))

        elif typ == '2d':
            # Create the coordinates
            x, y = np.meshgrid(x,y)
            x = x.T.flatten()
            y = y.T.flatten()
            z = z.flatten()

            ncols = (order + 1)**2
            G = np.zeros((z.size, ncols))
            ij = itertools.product(range(order+1), range(order+1)) # simply change here for uneven order?
            for k, (i,j) in enumerate(ij):
                G[:,k] = x**i * y**j
            w = np.sqrt(np.diag(weights)) # TODO add weights
            m, _, _, _ = np.linalg.lstsq(G, z)
            m = m.reshape(order+1,order+1)
            print m

        else:
            raise Exception('Typ must be 1d or 2d.')

        return m

    def  polyval(x=None, y=None, m=None, typ='1d'):
        """Values to two-dimensional polynomial fit. Based uppon code 
            provided by Joe Kington.
        """
        if typ == '1d':
            poly = np.poly1d(m)
            return poly(x)

        elif typ == '2d': #TODO: use polyval2d?
            x, y = np.meshgrid(x,y)
            x = x.T.flatten()
            y = y.T.flatten()
            return np.polynomial.polynomial.polyval2d(x, y, m)

#            order = int(np.sqrt(len(m))) - 1
#            ij = itertools.product(range(order+1), range(order+1))
#            z = np.zeros_like(x)
#            for a, (i,j) in zip(m, ij):
#                z += a * x**i * y**j
#            return z

        else:
            raise Exception('Typ must be 1d or 2d.')
 

    def clean_noisy(data, times, window, max_rms):
        """
        calculate a running RMS and remove noisy data

        window = in timestamps, sliding window dimension
        max_rms = flag points in a region with rms larger than max_rms times the rms of rmses
        
        return: an array of data dimensions with flags
        """
        # Check if data is empty
        if not data: 
            return []
        # loop over solution times
        rmses = np.zeros(shape=data.shape, dtype=np.float)
        for i, time in enumerate(times):
    
            # get data to smooth (values inside the time window)
            data_array = data[ np.where( abs(times - time) <= window / 2. ) ]
            rmses[i] = np.std(data_array)

        rms =  1.4826 * np.median( abs(rmses) )
        flags = abs(rmses) > max_rms * rms
 
        return flags
    
    def outlier_rej(vals, weights, axes, order=5, max_ncycles = 3, max_rms = 3., max_rms_noise = 2., replace = False):
        """
        Reject outliers using a running median
        val = the array (avg must be 0)
        weights = the weights to convert into flags
        axes = array with axes values (1d or 2d)
        max_ncycles = maximum number of cycles
        max_rms = number of rms times for outlier flagging
        replace = instead of flag it, replace the data point with the smoothed one
    
        return: flags array and final rms
        """
      
        if len(axes) == 1:
            axes[0] -= axes[0][0]
            typ = '1d'

        elif len(axes) == 2:
            axes[0] -= axes[0][0]
            axes[1] -= axes[1][0]
            typ = '2d'
        else:
            logging.error('FLAG operation can flag only along 1 or 2 axes. Given axes: '+str(axesToFlag))
            return

        if replace:
            orig_weights = np.copy(weights)
    
        for i in xrange(max_ncycles):

            # get polynomia and values
            if typ == '1d': 
                fit_sol = polyfit(axes[0], z=vals, w=weights, order=order, typ='1d')
                vals_detrend = vals - polyval(axes[0], m=fit_sol, typ='1d')

            elif typ == '2d': 
                fit_sol = polyfit(axes[0], axes[1], z=vals, w=weights, order=order, typ='2d')
                vals_detrend = vals - polyval(axes[0], axes[1], m=fit_sol, typ='2d').reshape( len( axes[0]), len(axes[1]) )
            print 'vals', vals
            print 'fit', polyval(axes[0], axes[1], m=fit_sol, typ='2d').reshape( len( axes[0]), len(axes[1]) )
            print 'detrend', vals_detrend

            # remove noisy regions of data
            #flag_noisy = clean_noisy(vals_detrend, time[ s ], window, max_rms_noise)
            #print 'noise flagging', float(sum(flag_noisy))/len(flag_noisy)
            #flags[ s ] = flag_noisy # add flags g (shape s=True) to global flags
            #s[ s ] = ~flag_noisy # new refined selection
            #vals_detrend = vals[ s ] - vals_smoothed[ ~flag_noisy ] # keep only vals satisfying s and g

            # all is flagged? break
            if (weights == 0).all():
                rms == 0.
                break

            # median calc
            rms =  1.4826 * np.median( np.abs(vals_detrend) )
            print 'rms', rms
    
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
        #plt.plot(axes[0], polyval(axes[0], m=fit_sol, typ='1d'), 'r-')
        #plt.savefig('test.png')

        # plot 2d
        import matplotlib.pyplot as plt
        plt.imshow(vals, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=5)
        plt.colorbar()
        plt.savefig('test2d.png')
        plt.clf()
        plt.imshow(polyval(axes[0], axes[1], m=fit_sol, typ='2d').reshape( len( axes[0]), len(axes[1]) ), origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=5)
        plt.colorbar()
        plt.savefig('test2d-smooth.png')
        plt.clf()
        plt.imshow(vals_detrend, origin='lower', interpolation="none", cmap=plt.cm.rainbow, aspect=5)
        plt.colorbar()
        plt.savefig('test2d-detrend.png')
        sys.exit(1)

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
        weights_re, re, rms_re = outlier_rej(re, weights, flagCoord, order, maxCycles, maxRms, maxRmsNoise, replace)
        weights_im, im, rms_im = outlier_rej(im, weights, flagCoord, order, maxCycles, maxRms, maxRmsNoise, replace)
        vals = np.arctan2(im, re)
        np.putmask(weights, weights_re == 0, 0)
        np.putmask(weights, weights_im == 0, 0)
        rms = np.sqrt(rms_re**2 + rms_im**2)

    elif solType == 'amplitude':
        weights, vals, rms = outlier_rej(np.log10(vals), weights, flagCoord, order, maxCycles, maxRms, maxRmsNoise, replace)
        vals == 10**vals

    else:
        weights, vals, rms = outlier_rej(vals, weights, flagCoord, order, maxCycles, maxRms, maxRmsNoise, replace)
    
    if percentFlagged(weights) == initPercentFlag:
        logging.debug('Percentage of data flagged/replaced (%s): None' % (removeKeys(coord, axisToFlag)))
    else: 
        logging.debug('Percentage of data flagged/replaced (%s): %.3f -> %.3f %% (rms: %.5f)' \
            % (removeKeys(coord, axisToFlag), initPercentFlag, percentFlagged(weights), rms))
    sys.exit(1)

#    outQueue.put([vals, flags, selection])
    return vals, weights, selection
        
            
def run( step, parset, H ):

    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    axesToFlag = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), 'time' )
    maxCycles = parset.getInt('.'.join(["LoSoTo.Steps", step, "MaxCycles"]), 5 )
    maxRms = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRms"]), 5. )
    maxRmsNoise = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRmsNoise"]), 5. )
    order = parset.getInt('.'.join(["LoSoTo.Steps", step, "Order"]), 3 )
    replace = parset.getBool('.'.join(["LoSoTo.Steps", step, "Replace"]), False )
    preflagzeros = parset.getBool('.'.join(["LoSoTo.Steps", step, "PreFlagZeros"]), False )
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 1 )
    
    if axesToFlag == []:
        logging.error("Please specify axis to flag. It must be a single one.")
        return 1

    # start processes for multi-thread
#    mpm = multiprocManager(ncpu, flag)

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Flagging soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        #sw = solWriter(soltab, useCache=True) # remember to flush!
        sw = solWriter(soltab) # remember to flush!

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
        axesToFlag = [axisToFlag for (coord, axisToFlag) in sorted(zip(sf.getAxesNames(),axesToFlag))]

        solType = sf.getType()

        # fill the queue (note that sf and sw cannot be put into a queue since they have file references)
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=axesToFlag, weight=True):
            #mpm.put([vals, weights, coord, solType, order, preflagzeros, maxCycles, maxRms, maxRmsNoise, replace, axesToFlag, selection])
            v, w, sel = flag(vals, weights, coord, solType, order, preflagzeros, maxCycles, maxRms, maxRmsNoise, replace, axesToFlag, selection)

#        mpm.wait()
        
#        for v, w, sel in mpm.get():
            sw.selection = sel
#            if replace:
#                # rewrite solutions (flagged values are overwritten)
#                sw.setValues(v, weight=False)
#            else:
#                # convert boolean flag to 01 float array (0->flagged)
#                # TODO: in this operation weight != 0,1 are lost
#                sw.setValues(w.astype(float), weight=True)
        
        #sw.flush()
        sw.addHistory('FLAG (over %s with %s sigma cut)' % (axesToFlag, maxRms))

        del sw
        del sf
        del soltab

    return 0
