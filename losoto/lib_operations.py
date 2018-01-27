# -*- coding: utf-8 -*-

# Some utilities for operations

import sys, math
import logging
from losoto.h5parm import h5parm
import multiprocessing
import numpy as np

class multiprocManager(object):

    class multiThread(multiprocessing.Process):
        """
        This class is a working thread which load parameters from a queue and
        return in the output queue
        """

        def __init__(self, inQueue, outQueue, funct):
            multiprocessing.Process.__init__(self)
            self.inQueue = inQueue
            self.outQueue = outQueue
            self.funct = funct

        def run(self):

            while True:
                parms = self.inQueue.get()

                # poison pill
                if parms is None:
                    self.inQueue.task_done()
                    break

                self.funct(*parms, outQueue=self.outQueue)
                self.inQueue.task_done()


    def __init__(self, procs=1, funct=None):
        """
        Manager for multiprocessing
        procs: number of processors
        funct: function to parallelize / note that the last parameter of this function must be the outQueue
        and it will be linked to the output queue
        """
        self.procs = procs
        self._threads = []
        self.inQueue = multiprocessing.JoinableQueue()
        self.outQueue = multiprocessing.Queue()
        self.runs = 0
        
        logging.debug('Spawning %i threads...' % self.procs)
        for proc in xrange(self.procs):
            t = self.multiThread(self.inQueue, self.outQueue, funct)
            self._threads.append(t)
            t.start()

    def put(self, args):
        """
        Parameters to give to the next jobs sent into queue
        """
        self.inQueue.put(args)
        self.runs += 1

    def get(self):
        """
        Return all the results as an iterator
        """
        # NOTE: do not use queue.empty() check which is unreliable
        # https://docs.python.org/2/library/multiprocessing.html
        for run in xrange(self.runs):
            yield self.outQueue.get()

    def wait(self):
        """
        Send poison pills to jobs and wait for them to finish
        The join() should kill all the processes
        """
        for t in self._threads:
            self.inQueue.put(None)

        # wait for all jobs to finish
        self.inQueue.join()


def reorderAxes( a, oldAxes, newAxes ):
    """
    Reorder axis of an array to match a new name pattern.

    Parameters
    ----------
    a : np array
        The array to transpose.
    oldAxes : list of str
        A list like ['time','freq','pol'].
        It can contain more axes than the new list, those are ignored.
        This is to pass to oldAxis the soltab.getAxesNames() directly even on an array from getValuesIter()
    newAxes : list of str
        A list like ['time','pol','freq'].

    Returns
    -------
    np array
        With axis transposed to match the newAxes list.
    """
    oldAxes = [ax for ax in oldAxes if ax in newAxes]
    idx = [ oldAxes.index(ax) for ax in newAxes ]
    return np.transpose(a, idx)


def removeKeys( dic, keys = [] ):
    """
    Remove a list of keys from a dict and return a new one.
    
    Parameters
    ----------
    dic : dcit
        The input dictionary.
    keys : list of str
        A list of arguments to remove or a string for single argument.

    Returns
    -------
    dict
        Dictionary with removed keys.
    """
    dicCopy = dict(dic)
    if type(keys) is str: keys = [keys]
    for key in keys:
        del dicCopy[key]
    return dicCopy


def normalize_phase(phase):
    """
    Normalize phase to the range [-pi, pi].
    
    Parameters
    ----------
    phase : array of float
        Phase to normalize.
    
    Returns
    -------
    array of float
        Normalized phases.
    """

    # Convert to range [-2*pi, 2*pi].
    out = np.fmod(phase, 2.0 * np.pi)
    # Remove nans
    nans = np.isnan(out)
    np.putmask(out, nans, 0)
    # Convert to range [-pi, pi]
    out[out < -np.pi] += 2.0 * np.pi
    out[out > np.pi] -= 2.0 * np.pi
    # Put nans back
    np.putmask(out, nans, np.nan)
    return out
