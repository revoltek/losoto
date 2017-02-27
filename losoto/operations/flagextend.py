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
        def extendFlag(flags, percent):
            #flags = flags.astype(np.int)
            if float(np.sum( flags ))/len(flags) > percent/100.:
                return 1
            else:
                return 0

        import scipy.ndimage
        initialPercent = 100.*(np.size(weights)-np.count_nonzero(weights))/np.size(weights)

        # if size=0 then extend to all 2*axis, this otherwise create issues with mirroring
        for i, s in enumerate(size):
            if s == 0: size[i] = 2*weights.shape[i]

        for cycle in xrange(cycles):
            flag = scipy.ndimage.filters.generic_filter((weights==0), extendFlag, size=size, mode='mirror', cval=0.0, origin=0, extra_keywords={'percent':percent})
            weights[ ( flag == 1 ) ] = 0
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
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 0 )
    if ncpu == 0:
        import multiprocessing
        ncpu = multiprocessing.cpu_count()
    
    if axesToExt == []:
        logging.error("Please specify at least one axis to extend flag.")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        # start processes for multi-thread
        mpm = multiprocManager(ncpu, flag)

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
                mpm.wait()
                return 1

        # fill the queue (note that sf and sw cannot be put into a queue since they have file references)
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=axesToExt, weight=True):
            mpm.put([weights, coord, axesToExt, selection, percent, size, cycles])

        mpm.wait()

        logging.info('Writing solutions')
        for w,sel in mpm.get():
            sw.selection = sel
            sw.setValues(w, weight=True) # convert back to np.float16

        sw.addHistory('FLAG EXTENDED (over %s)' % (str(axesToExt)))
        del sf
        del sw
    return 0
