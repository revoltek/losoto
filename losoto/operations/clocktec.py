#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an example operation for LoSoTo

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading CLOCKTEC module.')

def _run_parser(soltab, parser, step):
    clocksoltabOut = parser.getstr( step, 'clocksoltabOut', 'clock000' )
    tecsoltabOut = parser.getstr( step, 'tecsoltabOut', 'tec000' )
    offsetsoltabOut = parser.getstr( step, 'offsetsoltabOut', 'phase_offset000' )
    tec3rdsoltabOut = parser.getstr( step, 'tec3rdsoltabOut', 'tec3rd000' )
    flagBadChannels = parser.getbool( step, 'flagBadChannels', True )
    flagCut = parser.getfloat( step, 'flagCut', 5. )
    chi2cut = parser.getfloat( step, 'chi2cut', 3000. )
    combinePol = parser.getbool( step, 'combinePol', False )
    removePhaseWraps = parser.getbool( step, 'removePhaseWraps', True )
    fit3rdorder = parser.getbool( step, 'fit3rdorder', False )
    circular = parser.getbool( step, 'circular', False )
    reverse = parser.getbool( step, 'reverse', False )
    invertOffset = parser.getbool( step, 'invertOffset', False )
    nproc = parser.getint( step, 'nproc', 10 )

    parser.checkSpelling( step, soltab, ['tecsoltabOut', 'clocksoltabOut', 'offsetsoltabOut', 'tec3rdsoltabOut', 'flagBadChannels', 'flagCut', 'chi2cut', 'combinePol', 'removePhaseWraps', 'fit3rdorder', 'circular', 'reverse', 'invertOffset', 'nproc'])
    return run(soltab, tecsoltabOut, clocksoltabOut, offsetsoltabOut, tec3rdsoltabOut, flagBadChannels, flagCut, chi2cut, combinePol, removePhaseWraps, fit3rdorder, circular, reverse, invertOffset, nproc)


def run( soltab, tecsoltabOut='tec000', clocksoltabOut='clock000', offsetsoltabOut='phase_offset000', tec3rdsoltabOut='tec3rd000', flagBadChannels=True, flagCut=5., chi2cut=3000., combinePol=False, removePhaseWraps=True, fit3rdorder=False, circular=False, reverse=False, invertOffset=False, nproc=10 ):
    """
    Separate phase solutions into Clock and TEC.
    The Clock and TEC values are stored in the specified output soltab with type 'clock', 'tec', 'tec3rd'.

    Parameters
    ----------
    flagBadChannels : bool, optional
        Detect and remove bad channel before fitting, by default True.

    flagCut : float, optional


    chi2cut : float, optional


    combinePol : bool, optional
        Find a combined polarization solution, by default False.

    removePhaseWraps : bool, optional
        Detect and remove phase wraps, by default True.

    fit3rdorder : bool, optional
        Fit a 3rd order ionospheric ocmponent (usefult <40 MHz). By default False.

    circular : bool, optional
        Assume circular polarization with FR not removed. By default False.

    reverse : bool, optional
        Reverse the time axis. By default False.

    invertOffset : bool, optional
        Invert (reverse the sign of) the phase offsets. By default False. Set to True
        if you want to use them with the residuals operation.
    """
    import numpy as np
    from ._fitClockTEC import doFit

    logging.info("Clock/TEC separation on soltab: "+soltab.name)

    # some checks
    solType = soltab.getType()
    if solType != 'phase':
       logging.warning("Soltab type of "+soltab.name+" is: "+solType+" should be phase. Ignoring.")
       return 1

    # Collect station properties
    solset = soltab.getSolset()
    station_dict = solset.getAnt()
    stations = soltab.getAxisValues('ant')
    station_positions = np.zeros((len(stations), 3), dtype=np.float)
    for i, station_name in enumerate(stations):
        station_positions[i, 0] = station_dict[station_name][0]
        station_positions[i, 1] = station_dict[station_name][1]
        station_positions[i, 2] = station_dict[station_name][2]

    returnAxes=['ant','freq','pol','time']
    for vals, flags, coord, selection in soltab.getValuesIter(returnAxes=returnAxes,weight=True):

        if len(coord['ant']) < 2:
            logging.error('Clock/TEC separation needs at least 2 antennas selected.')
            return 1
        if len(coord['freq']) < 10:
            logging.error('Clock/TEC separation needs at least 10 frequency channels, preferably distributed over a wide range')
            return 1

        freqs=coord['freq']
        stations=coord['ant']
        times=coord['time']

        # get axes index
        axes=[i for i in soltab.getAxesNames() if i in returnAxes]

        # reverse time axes
        if reverse:
            vals = np.swapaxes(np.swapaxes(vals, 0, axes.index('time'))[::-1], 0, axes.index('time'))
            flags = np.swapaxes(np.swapaxes(flags, 0, axes.index('time'))[::-1], 0, axes.index('time'))

        result=doFit(vals,flags==0,freqs,stations,station_positions,axes,\
                         flagBadChannels=flagBadChannels,flagcut=flagCut,chi2cut=chi2cut,combine_pol=combinePol,removePhaseWraps=removePhaseWraps,fit3rdorder=fit3rdorder,circular=circular,n_proc=nproc)
        if fit3rdorder:
            clock,tec,offset,tec3rd=result
            if reverse:
                clock = clock[::-1,:]
                tec = tec[::-1,:]
                tec3rd = tec3rd[::-1,:]
        else:
            clock,tec,offset=result
            if reverse:
                clock = clock[::-1,:]
                tec = tec[::-1,:]
        if invertOffset:
            offset *= -1.0

        weights=tec>-5
        tec[np.logical_not(weights)]=0
        clock[np.logical_not(weights)]=0
        weights=np.float16(weights)

        if combinePol or not 'pol' in soltab.getAxesNames():
            tf_st = solset.makeSoltab('tec', soltabName = tecsoltabOut,
                             axesNames=['time', 'ant'], axesVals=[times, stations],
                             vals=tec[:,:,0],
                             weights=weights[:,:,0])
            tf_st.addHistory('CREATE (by CLOCKTECFIT operation)')
            tf_st = solset.makeSoltab('clock', soltabName = clocksoltabOut,
                             axesNames=['time', 'ant'], axesVals=[times, stations],
                             vals=clock[:,:,0]*1e-9,
                             weights=weights[:,:,0])
            tf_st.addHistory('CREATE (by CLOCKTECFIT operation)')
            tf_st = solset.makeSoltab('phase', soltabName = offsetsoltabOut,
                             axesNames=['ant'], axesVals=[stations],
                             vals=offset[:,0],
                             weights=np.ones_like(offset[:,0],dtype=np.float16))
            tf_st.addHistory('CREATE (by CLOCKTECFIT operation)')
            if fit3rdorder:
                tf_st = solset.makeSoltab('tec3rd', soltabName = tec3rdsoltabOut,
                                     axesNames=['time', 'ant'], axesVals=[times, stations],
                                     vals=tec3rd[:,:,0],
                                     weights=weights[:,:,0])
        else:
            tf_st = solset.makeSoltab('tec', soltabName = tecsoltabOut,
                             axesNames=['time', 'ant','pol'], axesVals=[times, stations, ['XX','YY']],
                             vals=tec,
                             weights=weights)
            tf_st.addHistory('CREATE (by CLOCKTECFIT operation)')
            tf_st = solset.makeSoltab('clock', soltabName = clocksoltabOut,
                             axesNames=['time', 'ant','pol'], axesVals=[times, stations, ['XX','YY']],
                             vals=clock*1e-9,
                             weights=weights)
            tf_st.addHistory('CREATE (by CLOCKTECFIT operation)')
            tf_st = solset.makeSoltab('phase', soltabName = offsetsoltabOut,
                             axesNames=['ant','pol'], axesVals=[stations, ['XX','YY']],
                             vals=offset,
                             weights=np.ones_like(offset,dtype=np.float16))
            tf_st.addHistory('CREATE (by CLOCKTECFIT operation)')
            if fit3rdorder:
                tf_st = solset.makeSoltab('tec3rd', soltabName = tec3rdsoltabOut,
                                     axesNames=['time', 'ant','pol'], axesVals=[times, stations, ['XX','YY']],
                                     vals=tec3rd,
                                     weights=weights)
    return 0
