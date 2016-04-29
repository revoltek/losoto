#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a extend flag procedure
# It can work in multi dimensional space and for each datum check if the surrounding data are flagged to a certain %, then flag also that datum
# The size of the surrounding footprint can be tuned
# WEIGHT: compliant

import logging
from losoto.operations_lib import *
import numpy as np

logging.debug('Loading FLAGEXTEND module.')

def flag(weights, coord, axesToExt, selection, percent=90, size=[0], cycles=3, outQueue=None):
        """
        Flag data if surreounded by other flagged data
        weights = the weights to convert into flags
        percent = percent of surrounding flagged point to extend the flag
        
        return: flags array and final rms
        """
        def extendFlag(weights, percent):
            if float(sum(weights))/len(weights) < 1.-(percent/100.):
                #print weights, float(sum(weights))/len(weights)
                return 1
            else:
                return 0

        import scipy.ndimage
        initialPercent = 100.*(np.size(weights)-np.count_nonzero(weights))/np.size(weights)

        # if size=0 then extend to all axis
        print size
        for i, s in enumerate(size):
            if s == 0: size[i] = weights.shape[i]
        print size

        for cycle in xrange(cycles):
            flag = scipy.ndimage.filters.generic_filter(weights, extendFlag, size=size, mode='mirror', cval=0.0, origin=0, extra_keywords={'percent':percent})
            weights[ np.where( flag == 1 ) ] = 0
            # no new flags
            if cycle != 0 and np.count_nonzero(flag) == oldFlagCount: break
            oldFlagCount = np.count_nonzero(flag)

        logging.debug('Percentage of data flagged (%s): %.3f -> %.3f %%' \
            % (removeKeys(coord, axesToExt), initialPercent, 100.*(np.size(weights)-np.count_nonzero(weights))/np.size(weights)))

        outQueue.put([weights, selection])
        
            
def run( step, parset, H ):

    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    axesToExt = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), ['freq','time'] )
    size = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "Size"]), [11,11] )
    percent = parset.getFloat('.'.join(["LoSoTo.Steps", step, "Percent"]), 50 )
    cycles = parset.getInt('.'.join(["LoSoTo.Steps", step, "Cycles"]), 3 )
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 1 )
    
    if axesToExt == []:
        logging.error("Please specify at least one axis to extend flag.")
        return 1

    # start processes for multi-thread
    mpm = multiprocManager(ncpu, flag)

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Extending flag on soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        for axisToExt in axesToExt:
            if axisToExt not in sf.getAxesNames():
                logging.error('Axis \"'+axisToExt+'\" not found.')
                return 1

        # fill the queue (note that sf and sw cannot be put into a queue since they have file references)
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=axesToExt, weight=True):
            # convert to float64 or numpy.ndimage complains
            mpm.put([weights.astype(np.float64), coord, axesToExt, selection, percent, size, cycles])

        mpm.wait()

        logging.info('Writing solutions')
        for w,sel in mpm.get():
            sw.selection = sel
            sw.setValues(w.astype(np.float16), weight=True) # convert back to np.float16

        sw.addHistory('FLAG EXTENDED (over %s)' % (str(axesToExt)))
        del sf
        del sw
    return 0
