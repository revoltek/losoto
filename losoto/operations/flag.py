#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a flagging procedure
# WEIGHT: flag-only compliant

import logging
from losoto.operations_lib import *
import numpy as np

logging.debug('Loading FLAG module.')

import multiprocessing
inQueue = multiprocessing.JoinableQueue()
outQueue = multiprocessing.Queue()

class multiThread(multiprocessing.Process):
    """
    This class is a working thread which load parameters from a queue and
    run the flagging on a chunk of data
    """

    def __init__(self, inQueue, outQueue):
        multiprocessing.Process.__init__(self)
        self.inQueue = inQueue
        self.outQueue = outQueue

    def run(self):

        while True:
            parms = self.inQueue.get()

            # poison pill
            if parms is None:
                self.inQueue.task_done()
                break

            self.flag(*parms)
            self.inQueue.task_done()

    def flag(self, vals, weights, coord, solType, preflagzeros, maxCycles, maxRms, maxRmsNoise, window, order, maxGap, replace, axisToFlag, selection):

        def smooth(data, times, window = 60., order = 1, max_gap = 5.*60. ):
            """
            Remove a trend from the data
            window = in timestamps, sliding window dimension
            order = 0: remove avg, 1: remove linear, 2: remove cubic
            max_gap = maximum allawed gap
        
            return: detrendized data array
            """
        
            final_data = np.copy(data)
        
            # loop over solution times
            for i, time in enumerate(times):
        
                # get data to smooth (values inside the time window)
                data_array = data[ np.where( abs(times - time) <= window / 2. ) ]
                data_offsets = times[ np.where( abs(times - time) <= window / 2. ) ] - time
        
                # check and remove big gaps in data
                if len( data_offsets ) > 1:
                  ddata_offsets = data_offsets[ 1 : ] - data_offsets[ : -1 ]
                  sel = np.where( ddata_offsets > max_gap )[0]
                  if len( sel ) > 0 :
                    min_data_index = 0
                    max_data_index = len( data_offsets )
                    this_time_index = np.where( abs( data_offsets ) == abs( data_offsets ).min() )[0]
                    # find min and max good indexes for this window
                    for s in sel:
                      if ( s < this_time_index ):
                        min_data_index = s + 1
                      if ( s >= this_time_index ):
                        max_data_index = s + 1
                        break
                    # redefine data arrays
                    data_array = data_array[ min_data_index : max_data_index ]
                    data_offsets = data_offsets[ min_data_index : max_data_index ]

                # smooth
                if len( data_array ) > 1:
                  dim = min( len( data_array ) - 2, order )
                  if ( dim == 0 ):
                    smooth_data = np.median( data_array )
                  else:
                    P = np.zeros( ( len( data_offsets ), dim + 1 ), dtype = data_offsets.dtype )
                    P[ : , 0 ] = 1.
                    if ( dim >= 1 ):
                        P[ : , 1 ] = data_offsets
                    if ( dim >= 2 ):
                        P[ : , 2 ] = data_offsets**2
                    Pt = np.transpose( P )
                    smooth_data = np.dot( np.linalg.inv( np.dot( Pt, P ) ), np.dot( Pt, data_array ) )[0]
                  final_data[i] = smooth_data
        
            return final_data
        ######################################

        def clean_noisy(data, times, window, max_rms):
            """
            calculate a running RMS and remove noisy data

            window = in timestamps, sliding window dimension
            max_rms = flag points in a region with rms larger than max_rms times the rms of rmses
            
            return: an array of data dimensions with flags
            """
            if len(data) == 0: return []
            # loop over solution times
            rmses = np.zeros(shape=data.shape, dtype=np.float)
            for i, time in enumerate(times):
        
                # get data to smooth (values inside the time window)
                data_array = data[ np.where( abs(times - time) <= window / 2. ) ]
                rmses[i] = np.std(data_array)

            rms =  1.4826 * np.median( abs(rmses) )
            flags = abs(rmses) > max_rms * rms
 
            return flags
        
        def outlier_rej(vals, weights, time, max_ncycles = 10, max_rms = 3., max_rms_noise=2., window = 60., order = 1, max_gap = 5.*60., replace = False):
            """
            Reject outliers using a running median
            val = the array (avg must be 0)
            weights = the weights to convert into flags
            time = array of seconds
            max_ncycles = maximum number of cycles
            max_rms = number of rms times for outlier flagging
            window, order, max_gap = see "smooth"
            replace = instead of flag it, replace the data point with the smoothed one
        
            return: flags array and final rms
            """
        
            flags = np.zeros(shape=weights.shape, dtype=np.bool)
            orig_flags = np.zeros(shape=weights.shape, dtype=np.bool)
            orig_flags[np.where(weights == 0.)] = True # initialize orig_flags to weights
        
            for i in xrange(max_ncycles):
        
                # smoothing (input with no flags!)
                s = ~orig_flags & ~flags # selecting non-flagged data
                vals_smoothed = smooth(vals[ s ], time[ s ], window, order, max_gap)
                vals_detrend = vals[ s ] - vals_smoothed
                
                # remove noisy regions of data
                flag_noisy = clean_noisy(vals_detrend, time[ s ], window, max_rms_noise)
                #print 'noise flagging', float(sum(flag_noisy))/len(flag_noisy)
                flags[ s ] = flag_noisy # add flags g (shape s=True) to global flags
                s[ s ] = ~flag_noisy # new refined selection
                vals_detrend = vals[ s ] - vals_smoothed[ ~flag_noisy ] # keep only vals satisfying s and g

                # all is flagged? break
                if (flags == True).all():
                    rms == 0.
                    break

                # median calc
                rms =  1.4826 * np.median( abs(vals_detrend) )
        
                # rejection  
                flag_outlier = abs(vals_detrend) > max_rms * rms
                #print 'rms flagging', float(sum(flag_outlier))/len(flag_outlier)

                flags[ s ] = flag_outlier
        
                # all is flagged? break
                if (flags == True).all():
                    rms == 0.
                    break
        
                # no flags? break
                if sum(flag_outlier) == 0 and sum(flag_noisy) == 0:
                    break
        
                # replace (outlier) flagged values with smoothed ones
                if replace:
                    new_vals = vals[ s ]
                    new_vals[ flag_outlier ] = vals_smoothed[ flag_outlier ]
                    vals[ s ] = new_vals

            return flags | orig_flags, vals, rms
        ########################################

        # check if everything flagged
        if np.count_nonzero(weights) == 0:
            logging.debug('Percentage of data flagged/replaced (%s): already completely flagged' % (removeKeys(coord, axisToFlag)))
            self.outQueue.put([vals, np.ones(shape=weights.shape, dtype=np.bool), selection])
            return

        if preflagzeros:
            if solType == 'amplitude': weights[np.where(vals == 1)] = 0
            else: weights[np.where(vals == 0)] = 0

        # if phase, then convert to real/imag, run the flagger on those, and convert back to pahses
        # best way to avoid unwrapping
        if solType == 'phase' or solType == 'scalarphase' or solType == 'rotation':
            re = 1. * np.cos(vals)
            im = 1. * np.sin(vals)
            flags_re, re, rms_re = outlier_rej(re, weights, coord[axisToFlag], maxCycles, maxRms, maxRmsNoise, window, order, maxGap, replace)
            flags_im, im, rms_im = outlier_rej(im, weights, coord[axisToFlag], maxCycles, maxRms, maxRmsNoise, window, order, maxGap, replace)
            vals = np.arctan2(im, re)
            flags = flags_re | flags_im
            rms = np.sqrt(rms_re**2 + rms_im**2)
            #flags, vals, rms = outlier_rej(unwrap(vals), weights, coord[axisToFlag], maxCycles, maxRms, maxRmsNoise, window, order, maxGap, replace)
            #vals = (vals+np.pi) % (2*np.pi) - np.pi
        elif solType == 'amplitude':
            flags, vals, rms = outlier_rej(np.log10(vals), weights, coord[axisToFlag], maxCycles, maxRms, maxRmsNoise, window, order, maxGap, replace)
            vals == 10**vals
        else:
            flags, vals, rms = outlier_rej(vals, weights, coord[axisToFlag], maxCycles, maxRms, maxRmsNoise, window, order, maxGap, replace)
        
        if (len(weights)-np.count_nonzero(weights))/float(len(weights)) == sum(flags)/float(len(flags)):
            logging.debug('Percentage of data flagged/replaced (%s): None' % (removeKeys(coord, axisToFlag)))
        else: 
            logging.debug('Percentage of data flagged/replaced (%s): %.3f -> %.3f %% (rms: %.5f)' \
                % (removeKeys(coord, axisToFlag), 100.*(len(weights)-np.count_nonzero(weights))/len(weights), 100.*sum(flags)/len(flags), rms))

        self.outQueue.put([vals, flags, selection])
        
            
def run( step, parset, H ):

    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    axisToFlag = parset.getString('.'.join(["LoSoTo.Steps", step, "Axis"]), 'time' )
    maxCycles = parset.getInt('.'.join(["LoSoTo.Steps", step, "MaxCycles"]), 5 )
    maxRms = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRms"]), 5. )
    maxRmsNoise = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRmsNoise"]), 5. )
    window = parset.getFloat('.'.join(["LoSoTo.Steps", step, "Window"]), 100 )
    order = parset.getInt('.'.join(["LoSoTo.Steps", step, "Order"]), 1 )
    maxGap = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxGap"]), 5*60 )
    replace = parset.getBool('.'.join(["LoSoTo.Steps", step, "Replace"]), False )
    preflagzeros = parset.getBool('.'.join(["LoSoTo.Steps", step, "PreFlagZeros"]), False )
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 1 )
    
    if axisToFlag == '':
        logging.error("Please specify axis to flag. It must be a single one.")
        return 1

    if order > 2 or order < 0:
        logging.error("Order must be 0 (mean), 1 (linear), 2 (cubic)")
        return 1

    # start processes for multi-thread
    logging.debug('Spowning %i threads...' % ncpu)
    for i in range(ncpu):
        t = multiThread(inQueue, outQueue)
        t.start()

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Flagging soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab, useCache=True) # remember to flush!

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        if axisToFlag not in sf.getAxesNames():
            logging.error('Axis \"'+axis+'\" not found.')
            return 1

        solType = sf.getType()

        # fill the queue (note that sf and sw cannot be put into a queue since they have file references)
        runs = 0
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=axisToFlag, weight=True):
            runs += 1
            inQueue.put([vals, weights, coord, solType, preflagzeros, maxCycles, maxRms, maxRmsNoise, window, order, maxGap, replace, axisToFlag, selection])

        # add poison pills to kill processes
        for i in range(ncpu):
            inQueue.put(None)

        # wait for all jobs to finish
        inQueue.join()
        
        # writing back the solutions
        # NOTE: do not use queue.empty() check which is unreliable
        # https://docs.python.org/2/library/multiprocessing.html
        for i in range(runs):
            q = outQueue.get()
            v,f,sel = q
            sw.selection = sel
            if replace:
                # rewrite solutions (flagged values are overwritten)
                sw.setValues(v, weight=False)
            else:
                # convert boolean flag to 01 float array (0->flagged)
                # TODO: in this operation weight != 0,1 are lost
                sw.setValues((~f).astype(float), weight=True)

        sw.flush()

        sw.addHistory('FLAG (over %s with %s sigma cut)' % (axisToFlag, maxRms))
    return 0
