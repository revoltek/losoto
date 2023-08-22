#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging
from losoto.operations._faraday_timestep import _run_timestep

logging.debug('Loading FRjump module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'rotationmeasure002' )
    soltabPhase = parser.getstr( step, 'soltabPhase', 'phase000')
    clipping = np.array(parser.getarray(step, 'clipping', [0,1e9]),dtype=float)
    frequencies = np.array(parser.getarray(step, 'frequencies', []),dtype=float)

    parser.checkSpelling( step, soltab, ['soltabOut','soltabPhase','clipping','frequencies'])
    return run(soltab, soltabOut,clipping, soltabPhase,frequencies)

def costfunctionRM(RM, wav, phase):
    return np.sum(abs(np.cos(2.*RM[0]*wav*wav) - np.cos(phase)) + abs(np.sin(2.*RM[0]*wav*wav) - np.sin(phase)))

def getPhaseWrapBase(wavels):
    """
    freqs: frequency grid of the data
    return the step size from a local minima (2pi phase wrap) to the others [0]: TEC, [1]: clock
    """
    wavels = np.array(wavels)
    nF = wavels.shape[0]
    A = np.zeros((nF, 1), dtype=float)
    A[:, 0] = 2*wavels**2
    steps = np.dot(np.dot(np.linalg.inv(np.dot(A.T, A)), A.T), 2 * np.pi * np.ones((nF, ), dtype=float))
    return steps

def getPhaseWrapBase_TEC(freqs):
    """
    freqs: frequency grid of the data
    return the step size from a local minima (2pi phase wrap) to the others [0]: TEC, [1]: clock
    """
    freqs = np.array(freqs)
    nF = freqs.shape[0]
    A = np.zeros((nF, 2), dtype=float)
    A[:, 1] = freqs * 2 * np.pi * 1e-9
    A[:, 0] = -8.44797245e9 / freqs
    steps = np.dot(np.dot(np.linalg.inv(np.dot(A.T, A)), A.T), 4 * np.pi * np.ones((nF, ), dtype=float))
    return steps[0]

def dejump(vals,wavels,dotec=False):
    if dotec:
        c = 2.99792458e8
        freqs = c/wavels
        phasewrap = np.abs(getPhaseWrapBase_TEC(freqs))
    else:
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

def run( soltab, soltabOut,clipping,soltabPhase,frequencies):
    """
    EXPERIMENTAL
    'Dejumps' the Faraday solutions. Because Faraday rotation is a rotation, there are generally multiple possible values for the rotation measure
    that yield a similar rotation angle - but are offset from the main trend. This code will fix this.

    Parameters
    ----------
    
    soltabOut : str, optional
        output table name (same solset), by deault "rotationmeasure001".
        
    clipping : arr, optional
        Refers to the frequency range that is used in the Faraday step. This is important to get a good estimate for the difference in RM

    soltabPhase : str, optional
        name of soltab that contains the phases. The only reason we need this, is for the frequency axis (so it really doesnt matter what's in 
        this soltab. Default: phase000

    frequencies : arr, optional
        if the frequency axis is missing, you need to give the freq range used for TEC/FR fitting manually here.
    """

    import numpy as np
    import scipy.optimize

    c = 2.99792458e8

    vals = soltab.getValues(retAxesVals=False)
    ants = soltab.getAxisValues('ant')
    times = soltab.getAxisValues('time')
    solset = soltab.getSolset()
    if len(frequencies) == 0:
        phases = solset.getSoltab(soltabPhase)
        freqs = phases.getValues()[1]['freq']
    else:
        freqs = np.arange(*frequencies)
    selection = (freqs > clipping[0]) * (freqs < clipping[1]) # Only take the frequencies used for fitting
    selected_freqs = freqs[selection]

    wavels = c/selected_freqs # in meters

    if soltabOut not in solset.getSoltabNames():
        soltabout = solset.makeSoltab(soltab.getType(),soltabName=soltabOut, 
                            axesNames=['ant','time'], axesVals=[ants,times],
                            vals = np.zeros((len(ants),len(times))),
                            weights = np.ones((len(ants),len(times))))
    else:
        soltabout = solset.getSoltab(soltabOut)
    soltabout.addHistory('Created by the FRJUMP operation from %s.'%soltab.name)

    dotec = (soltab.getType()=='tec')

    # Iterate through all antennas
    
    for vals,weights,coord,_ in soltab.getValuesIter(returnAxes='time',weight=True):
        antname = coord['ant']
        newvals = dejump(vals,wavels,dotec)
        logging.debug(f'Doing antenna {antname}')
        soltabout.setSelection(ant=coord['ant'],time=coord['time'])          
        soltabout.setValues(np.expand_dims(newvals,axis=1))
    


    return 0
