#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a flagging procedure

import logging
from operations_lib import *

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
            smooth_data = np.dot( np.linalg.inv( np.dot( Pt, P ) ), np.dot( Pt, data_array ) )[ 0 ]
          n = np.where( times == time )[0]
          final_data[n] = smooth_data

    return final_data


def outlier_rej(val, time, max_ncycles = 10, max_rms = 3., window = 60., order = 1, max_gap = 5.*60.):
    """
    Reject outliers using a running median
    val = the array (avg must be 0)
    time = array of seconds
    max_ncycles = maximum number of cycles
    max_rms = number of rms times for outlier flagging
    window, order, max_gap = see "smooth"

    return: flags array and final rms
    """
    flags = np.zeros(shape=val.shape, dtype=np.bool)
    val_detrend = np.zeros(shape=val.shape)

    for i in xrange(max_ncycles):

        # smoothing (input with no flags!)
        val_smoothed = smooth(val[~flags], time[~flags], window, order, max_gap)
        val_detrend[~flags] = val[~flags] - val_smoothed
        
        # median calc
        rms =  1.4826 * np.median(abs(val_detrend[~flags]))

        # rejection  
        flags[ ~flags ] = np.logical_or(flags[~flags], abs(val_detrend[~flags]) > max_rms * rms )

        if (flags == True).all():
            rms == 0.
            break

        # median calc
        this_rms =  1.4826 * np.median(abs(val_detrend[~flags]))
        if rms - this_rms == 0.:
            break

    return flags, rms


def run( step, parset, H ):

    import numpy as np
    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    axisToFlag = parset.getString('.'.join(["LoSoTo.Steps", step, "Axis"]), '' )
    maxCycles = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxCycles"]), 5. )
    maxRms = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxRms"]), 5. )
    window = parset.getFloat('.'.join(["LoSoTo.Steps", step, "Window"]), 10. )
    order = parset.getInt('.'.join(["LoSoTo.Steps", step, "Order"]), 1 )
    maxGap = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxGap"]), 5.*60. )
    
    if axisToFlag == '':
        logging.error("Please specify axis to flag. It must be a single one.")
        return 1
    if order > 2 or order < 0:
        logging.error("Order must be 0 (mean), 1 (linear), 2 (cubic)")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Smoothing soltab: "+soltab.name)

        sf.setSelection(ant=ants, pol=pols, dir=dirs)

        if axisToFlag not in sf.getAxes():
            logging.error('Axis \"'+axis+'\" not found.')
            return 1

        for vals, coord in sf.getValuesIter(returnAxes=axisToFlag):

            # if phase, then unwrap and flag
            if sf.getType() == 'phase' or sf.getType() == 'scalarphase':
                vals_unwrap = unwrap_fft(vals)
                flags, rms = outlier_rej(vals_unwrap, coord[axisToFlag], maxCycles, maxRms, window, order, maxGap)
            else:
                flags, rms = outlier_rej(vals, coord[axisToFlag], maxCycles, maxRms, window, order, maxGap)

            logging.info('Final rms: '+str(rms))

            # writing back the solutions
            coord = removeKeys(coord, axisToFlag)
            sw.setSelection(**coord)
            # convert boolean flag to 01 binary array (0->flagged)
            sw.setValues((~flags).astype(int), weigth=True)

        sw.addHistory('FLAG (over %s with %s sigma cut)' % (axisToFlag, maxRms))
    return 0
