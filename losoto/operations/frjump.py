#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging
from losoto.operations._faraday_timestep import _run_timestep
import astropy.constants as c

logging.debug('Loading FRjump module.')


def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'rotationmeasure002' )
    soltabPhase = parser.getstr( step, 'soltabPhase', 'phase000')
    clipping = parser.getarray(step, 'clipping', [0,1e9])

    parser.checkSpelling( step, soltab, ['soltabOut'])
    return run(soltab, soltabOut,clipping, soltabPhase)

def costfunctionRM(RM, wav, phase):
    return np.sum(abs(np.cos(2.*RM[0]*wav*wav) - np.cos(phase)) + abs(np.sin(2.*RM[0]*wav*wav) - np.sin(phase)))

def getPhaseWrapBase(wavels):
    """
    freqs: frequency grid of the data
    return the step size from a local minima (2pi phase wrap) to the others [0]: TEC, [1]: clock
    """
    wavels = np.array(wavels)
    nF = wavels.shape[0]
    A = np.zeros((nF, 1), dtype=np.float)
    A[:, 0] = 2*wavels**2
    steps = np.dot(np.dot(np.linalg.inv(np.dot(A.T, A)), A.T), 2 * np.pi * np.ones((nF, ), dtype=np.float))
    return steps

def dejump(vals,wavels):
    phasewrap = getPhaseWrapBase(wavels)
    diffs = np.diff(vals)
    
    # Fix the jumps in diff-space
    for i in range(len(diffs)):
        while np.abs(diffs[i]) > phasewrap/2.:
            if diffs[i] > phasewrap/2.:
                diffs[i] -= phasewrap
            elif diffs[i] < -phasewrap/2.:
                diffs[i] += phasewrap
                
    newvals = np.zeros_like(vals)
    newvals[0] = vals[0]
    for i in range(len(diffs)):
        newvals[i+1] = newvals[i] + diffs[i] # Here, we go back to normal-space
    
    delt_median = lambda nv: np.median(nv) - np.median(vals) # make sure that the median is within 1 phasejump
    while np.abs(delt_median(newvals)) > phasewrap/2.:
        if delt_median(newvals) > phasewrap/2.:
            newvals -= phasewrap
        elif delt_median(newvals) < -phasewrap/2.:
            newvals += phasewrap
    return newvals

def run( soltab, soltabOut,clipping,soltabPhase):
    """
    Faraday rotation extraction from either a rotation table or a circular phase (of which the operation get the polarisation difference).

    Parameters
    ----------
    
    soltabOut : str, optional
        output table name (same solset), by deault "rotationmeasure000".
        
    refAnt : str, optional
        Reference antenna, by default the first.

    maxResidual : float, optional
        Max average residual in radians before flagging datapoint, by default 1. If 0: no check.

    """
    import numpy as np
    import scipy.optimize

    c = 2.99792458e8

    vals = soltab.getValues(retAxesVals=False)
    ants = soltab.getAxisValues('ant')
    times = soltab.getAxisValues('time')
    solset = soltab.getSolset()
    phases = solset.getSoltab(soltabPhase)
    freqs = phases.getValues()[1]['freq']
    selection = (freqs > clipping[0]) * (freqs < clipping[1]) # Only take the frequencies used for fitting
    selected_freqs = freqs[selection]

    wavels = c/selected_freqs # in meters

    if soltabOut not in solset.getSoltabNames():
        soltabout = solset.makeSoltab('rotationmeasure',soltabName=soltabOut, 
                            axesNames=['ant','time'], axesVals=[ants,times],
                            vals = np.zeros((len(ants),len(times))),
                            weights = np.ones((len(ants),len(times))))
    else:
        soltabout = solset.getSoltab(soltabOut)
    soltabout.addHistory('Created by the FRJUMP operation from %s.'%soltab.name)


    # Iterate through all antennas
    
    for vals,weights,coord,_ in soltab.getValuesIter(returnAxes='time',weight=True):
        antname = coord['ant']
        newvals = dejump(vals,wavels)
        logging.info(f'Doing antenna {antname}')
        soltabout.setSelection(ant=coord['ant'],time=coord['time'])          
        soltabout.setValues(np.expand_dims(newvals,axis=1))
    


    return 0
