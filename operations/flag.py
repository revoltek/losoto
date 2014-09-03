#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a flagging procedure

import logging
from operations_lib import *
import numpy as np

logging.debug('Loading FLAG module.')


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
    for time in times:

        # get data to smooth (values inside the time window)
        data_array = data[ np.where( abs(times - time) <= window / 2. ) ]
        data_offsets = times[ np.where( abs(times - time) <= window / 2. ) ] - time

        # check and remove big gaps in data
        if ( len( data_offsets ) > 1 ):
          ddata_offsets = data_offsets[ 1 : ] - data_offsets[ : -1 ]
          sel = np.where( ddata_offsets > max_gap )[0]
          if ( len( sel ) > 0 ):
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
          final_data[np.where( times == time )] = smooth_data

    return final_data


def outlier_rej(vals, weights, time, max_ncycles = 10, max_rms = 3., window = 60., order = 1, max_gap = 5.*60., replace = False):
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
    orig_flags[np.where(weights == 0.)] == True # initialize orig_flags to weights

    for i in xrange(max_ncycles):

        # smoothing (input with no flags!)
        s = ~orig_flags & ~flags # selecting non-flagged data
        vals_smoothed = smooth(vals[ s ], time[ s ], window, order, max_gap)
        vals_detrend = vals[ s ] - vals_smoothed
        
        # median calc
        rms =  1.4826 * np.median( abs(vals_detrend) )

        # rejection  
        new_flags = abs(vals_detrend) > max_rms * rms
        flags[ s ] = new_flags

        # all is flagged? break
        if (flags == True).all():
            rms == 0.
            break

        # median calc
        this_rms =  1.4826 * np.median( abs(vals_detrend[ ~new_flags ]) )

        # no flags? break
        if rms - this_rms == 0.:
            break

        # replace flagged values with smoothed ones
        if replace:
            new_vals = vals[ s ]
            new_vals[ new_flags ] = vals_smoothed[ new_flags ]
            vals[ s ] = new_vals

    return flags, vals, rms


def run( step, parset, H ):

    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    axisToFlag = parset.getString('.'.join(["LoSoTo.Steps", step, "Axis"]), '' )
    maxCycles = parset.getInt('.'.join(["LoSoTo.Steps", step, "MaxCycles"]), 5 )
    maxRms = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRms"]), 5. )
    window = parset.getInt('.'.join(["LoSoTo.Steps", step, "Window"]), 10 )
    order = parset.getInt('.'.join(["LoSoTo.Steps", step, "Order"]), 1 )
    maxGap = parset.getInt('.'.join(["LoSoTo.Steps", step, "MaxGap"]), 5*60 )
    replace = parset.getBool('.'.join(["LoSoTo.Steps", step, "Replace"]), False )
    
    if axisToFlag == '':
        logging.error("Please specify axis to flag. It must be a single one.")
        return 1

    if order > 2 or order < 0:
        logging.error("Order must be 0 (mean), 1 (linear), 2 (cubic)")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Flagging soltab: "+soltab._v_name)

        sf.setSelection(ant=ants, pol=pols, dir=dirs)

        if axisToFlag not in sf.getAxesNames():
            logging.error('Axis \"'+axis+'\" not found.')
            return 1

        for vals, weights, coord in sf.getValuesIter(returnAxes=axisToFlag, weight=True):

            # if phase, then unwrap, flag and wrap again
            if sf.getType() == 'phase' or sf.getType() == 'scalarphase' or sf.getType() == 'rotation':
                flags, vals, rms = outlier_rej(unwrap_fft(vals), weights, coord[axisToFlag], maxCycles, maxRms, window, order, maxGap, replace)
                vals = (vals+np.pi) % (2*np.pi) - np.pi
            else:
                flags, vals, rms = outlier_rej(vals, weights, coord[axisToFlag], maxCycles, maxRms, window, order, maxGap, replace)

            logging.debug('Percentage of data flagged/replaced (%s): %.3f -> %.3f %%' \
                    % (removeKeys(coord, axisToFlag), 100.*len(np.where(weights==0)[0])/len(weights), 100.*sum(flags)/len(flags)))

            # writing back the solutions
            coord = removeKeys(coord, axisToFlag)
            sw.setSelection(**coord)
            if replace:
                # rewrite solutions (flagged values are overwritten)
                sw.setValues(vals, weight=False)
            else:
                # convert boolean flag to 01 binary array (0->flagged)
                sw.setValues((~flags).astype(int), weight=True)

        sw.addHistory('FLAG (over %s with %s sigma cut)' % (axisToFlag, maxRms))
    return 0
