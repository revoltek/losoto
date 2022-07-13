from losoto.lib_operations import *
from losoto._logging import logger as logging
import multiprocessing as mp
import numpy as np
import scipy.optimize

def costfunctionRM(RM, wav, phase):
    return np.sum(abs(np.cos(2.*RM[0]*wav*wav) - np.cos(phase)) + abs(np.sin(2.*RM[0]*wav*wav) - np.sin(phase)))

def _run_timestep(t,coord_rr,coord_ll,weights,vals,solType,coord,maxResidual):
    c = 2.99792458e8
    if solType == 'phase':
        idx       = ((weights[coord_rr,:] != 0.) & (weights[coord_ll,:] != 0.))
        freq      = np.copy(coord['freq'])[idx]
        phase_rr  = vals[coord_rr,:][idx]
        phase_ll  = vals[coord_ll,:][idx]
        # RR-LL to be consistent with BBS/NDPPP
        phase_diff  = (phase_rr - phase_ll)      # not divide by 2 otherwise jump problem, then later fix this
    else: # rotation table
        idx        = ((weights[:] != 0.) & (weights[:] != 0.))
        freq       = np.copy(coord['freq'])[idx]
        phase_diff = 2.*vals[:][idx] # a rotation is between -pi and +pi

    if len(freq) < 20:
        fitresultrm_wav = [0]
        weight = 0
        logging.warning('No valid data found for Faraday fitting for antenna: '+coord['ant']+' at timestamp '+str(t))
    else:
        # if more than 1/4 of chans are flagged
        if (len(idx) - len(freq))/float(len(idx)) > 1/4.:
            logging.debug('High number of filtered out data points for the timeslot %i: %i/%i' % (t, len(idx) - len(freq), len(idx)) )

        wav = c/freq

        #fitresultrm_wav, success = scipy.optimize.leastsq(rmwavcomplex, [fitrmguess], args=(wav, phase_diff))
        ranges = slice(-0.1, 0.1, 2e-4)
        fitresultrm_wav = scipy.optimize.brute(costfunctionRM, (ranges,), finish=scipy.optimize.leastsq, args=(wav, phase_diff))        

        # fractional residual
        residual = np.nanmean(np.abs(np.mod((2.*fitresultrm_wav*wav*wav)-phase_diff + np.pi, 2.*np.pi) - np.pi))
        if maxResidual == 0 or residual < maxResidual:
            fitrmguess = fitresultrm_wav[0] # Don't know what this is for...
            weight = 1
        else:       
            # high residual, flag
            logging.warning('Bad solution for ant: '+coord['ant']+' (time: '+str(t)+', residual: '+str(residual)+').')
            weight = 0

    return fitresultrm_wav[0],weight
