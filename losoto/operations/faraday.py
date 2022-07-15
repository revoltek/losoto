#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging
import multiprocessing as mp
from losoto.operations._faraday_timestep import _run_timestep

logging.debug('Loading FARADAY module.')


def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'rotationmeasure000' )
    refAnt = parser.getstr( step, 'refAnt', '')
    maxResidual = parser.getfloat( step, 'maxResidual', 1. )
    ncpu = parser.getint( step, 'ncpu', 0)

    parser.checkSpelling( step, soltab, ['soltabOut', 'refAnt', 'maxResidual','ncpu'])
    return run(soltab, soltabOut, refAnt, maxResidual, ncpu)

def costfunctionRM(RM, wav, phase):
    return np.sum(abs(np.cos(2.*RM[0]*wav*wav) - np.cos(phase)) + abs(np.sin(2.*RM[0]*wav*wav) - np.sin(phase)))

def run( soltab, soltabOut='rotationmeasure000', refAnt='', maxResidual=1.,ncpu=0):
    logging.info(f'ncpus: {ncpu}')
    """
    Faraday rotation extraction from either a rotation table or a circular phase (of which the operation get the polarisation difference).

    Parameters
    ----------
    
    soltabOut : str, optional
        output table name (same solset), by deault "rotationmeasure000".
        
    refAnt : str, optional
        Reference antenna, by default 'auto'.

    maxResidual : float, optional
        Max average residual in radians before flagging datapoint, by default 1. If 0: no check.

    """
    import numpy as np
    import scipy.optimize

    rmwavcomplex = lambda RM, wav, y: abs(np.cos(2.*RM[0]*wav*wav)  - np.cos(y)) + abs(np.sin(2.*RM[0]*wav*wav)  - np.sin(y))
    c = 2.99792458e8

    logging.info("Find FR for soltab: "+soltab.name)

    # input check
    solType = soltab.getType()
    if solType == 'phase':
        returnAxes = ['pol','freq','time']
        if 'RR' in soltab.getAxisValues('pol') and 'LL' in soltab.getAxisValues('pol'):
            coord_rr = np.where(soltab.getAxisValues('pol') == 'RR')[0][0]
            coord_ll = np.where(soltab.getAxisValues('pol') == 'LL')[0][0]
        elif 'XX' in soltab.getAxisValues('pol') and 'YY' in soltab.getAxisValues('pol'):
            logging.warning('Linear polarization detected, LoSoTo assumes XX->RR and YY->LL.')
            coord_rr = np.where(soltab.getAxisValues('pol') == 'XX')[0][0]
            coord_ll = np.where(soltab.getAxisValues('pol') == 'YY')[0][0]
        else:
            logging.error("Cannot proceed with Faraday estimation with polarizations: "+str(coord['pol']))
            return 1
    elif solType == 'rotation':
        returnAxes = ['freq','time']
        coord_rr = None
        coord_ll = None
    else:
       logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase or rotation. Ignoring.")
       return 1

    if refAnt == '': refAnt = 'auto'
    elif refAnt != 'closest' and refAnt != 'auto' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.warning('Reference antenna '+refAnt+' not found. Using: atomatic search.')
        refAnt = 'auto'

    # times and ants needs to be complete or selection is much slower
    times = soltab.getAxisValues('time')
    ants = soltab.getAxisValues('ant')

    # create new table
    solset = soltab.getSolset()
    soltabout = solset.makeSoltab('rotationmeasure', soltabName = soltabOut,
                             axesNames=['ant','time'], axesVals=[ants, times],
                             vals=np.zeros((len(ants),len(times))),
                             weights=np.ones((len(ants),len(times))))
    soltabout.addHistory('Created by FARADAY operation from %s.' % soltab.name)

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=returnAxes, weight=True, refAnt=refAnt):

        if len(coord['freq']) < 10:
            logging.error('Faraday rotation estimation needs at least 10 frequency channels, preferably distributed over a wide range.')
            return 1

        # reorder axes
        vals = reorderAxes( vals, soltab.getAxesNames(), returnAxes )
        weights = reorderAxes( weights, soltab.getAxesNames(), returnAxes )
        weights[np.isnan(vals)] = 0.

        fitrm = np.zeros(len(times))
        fitweights = np.ones(len(times)) # all unflagged to start
        fitrmguess = 0.001 # good guess

        if not coord['ant'] == refAnt and not (refAnt == 'auto' and coord['ant'] == soltab.cacheAutoRefAnt):
            logging.debug('Working on ant: '+coord['ant']+'...')

            if (weights == 0.).all() == True:
                logging.warning('Skipping flagged antenna: '+coord['ant'])
                fitweights[:] = 0
            else:
                if solType == 'phase':
                    weightsliced = [weights[:,:,t] for t,_ in enumerate(times)]
                    valsliced = [vals[:,:,t] for t,_ in enumerate(times)]
                else: # rotation table
                    weightsliced = [weights[:,t] for t,_ in enumerate(times)]
                    valsliced = [vals[:,t] for t,_ in enumerate(times)]

                tuples = [(t,coord_rr,coord_ll,wt,vl,solType,coord,maxResidual) for t,wt,vl in zip(list(np.arange(len(times))), weightsliced, valsliced)]
                if ncpu == 0:
                    ncpu = mp.cpu_count()
                with mp.Pool(ncpu) as pool:
                    fitrm,fitweights = zip(*pool.starmap(_run_timestep,tuples))

        soltabout.setSelection(ant=coord['ant'], time=coord['time'])
        soltabout.setValues( np.expand_dims(fitrm, axis=1) )
        soltabout.setValues( np.expand_dims(fitweights, axis=1), weight=True )

    return 0
