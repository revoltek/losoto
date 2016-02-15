#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a extend flag procedure
# WEIGHT: compliant

import logging
from losoto.operations_lib import *
import numpy as np

logging.debug('Loading FLAGEXTEND module.')

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

    def flag(self, weights, coord, axesToExt, selection, percent=90, size=11, cycles=3):
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
        for cycle in xrange(cycles):
            flag = scipy.ndimage.filters.generic_filter(weights, extendFlag, size=size, mode='mirror', cval=0.0, origin=0, extra_keywords={'percent':percent})
            weights[ np.where( flag == 1 ) ] = 0
            # no new flags
            if cycle != 0 and np.count_nonzero(flag) == oldFlagCount: break
            oldFlagCount = np.count_nonzero(flag)

        logging.debug('Percentage of data flagged (%s): %.3f -> %.3f %%' \
            % (removeKeys(coord, axesToExt), initialPercent, 100.*(np.size(weights)-np.count_nonzero(weights))/np.size(weights)))

        self.outQueue.put([weights, selection])
        
            
def run( step, parset, H ):

    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    axesToExt = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), ['freq','time'] )
    size = parset.getInt('.'.join(["LoSoTo.Steps", step, "Size"]), 11 )
    percent = parset.getFloat('.'.join(["LoSoTo.Steps", step, "Percent"]), 50 )
    cycles = parset.getInt('.'.join(["LoSoTo.Steps", step, "Cycles"]), 3 )
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 1 )
    
    if axesToExt == []:
        logging.error("Please specify at least one axis to extend flag.")
        return 1

    # start processes for multi-thread
    logging.debug('Spowning %i threads...' % ncpu)
    for i in xrange(ncpu):
        t = multiThread(inQueue, outQueue)
        t.start()

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
        runs = 0
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=axesToExt, weight=True):
            runs += 1
            # convert to float64 or numpy.ndimage complains
            inQueue.put([weights.astype(np.float64), coord, axesToExt, selection, percent, size, cycles])

        # add poison pills to kill processes
        for i in xrange(ncpu):
            inQueue.put(None)

        # wait for all jobs to finish
        inQueue.join()
        
        # writing back the solutions
        # NOTE: do not use queue.empty() check which is unreliable
        # https://docs.python.org/2/library/multiprocessing.html
        logging.info('Writing solutions')
        for i in xrange(runs):
            q = outQueue.get()
            w,sel = q
            sw.selection = sel
            sw.setValues(w.astype(np.float16), weight=True) # convert back to np.float16

        sw.addHistory('FLAG EXTENDED (over %s)' % (str(axesToExt)))
        del sf
        del sw
    return 0
