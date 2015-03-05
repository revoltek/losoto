#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an example operation for LoSoTo

import logging
from operations_lib import *
from fitClockTEC import doFit
logging.debug('Loading CLOCKTEC module.')

def run( step, parset, H ):
    """
    Separate phase solutions into Clock and TEC.

    Phase solutions are assumed to be stored in solsets of the H5parm file, one
    solset per field.

    The Clock and TEC values are stored in the specified output soltab with type 'clock' and 'tec'.

    """
    import numpy as np
    from h5parm import solFetcher, solWriter

    # get involved solsets using local step values or global values or all
    soltabs = getParSoltabs( step, parset, H )

    flagBadChannels = parset.getBool('.'.join(["LoSoTo.Steps", step, "FlagBadChannels"]), True )
    flagCut = parset.getFloat('.'.join(["LoSoTo.Steps", step, "FlagCut"]), 1.5 )
    chi2cut = parset.getFloat('.'.join(["LoSoTo.Steps", step, "Chi2cut"]), 30000. )
    combinePol = parset.getBool('.'.join(["LoSoTo.Steps", step, "CombinePol"]), False )
    fitOffset = parset.getBool('.'.join(["LoSoTo.Steps", step, "FitOffset"]), False )
    removePhaseWraps=parset.getBool('.'.join(["LoSoTo.Steps", step, "RemovePhaseWraps"]), True )
    # do something on every soltab (use the openSoltab LoSoTo function)
    #for soltab in openSoltabs( H, soltabs ):
    for soltabname in soltabs:
        solsetname=soltabname.split('/')[0]
        soltab=H.getSoltab(solset=solsetname, soltab=soltabname.split('/')[1])
        logging.info("--> Working on soltab: "+soltab._v_name)
        t = solFetcher(soltab)
        tw = solWriter(soltab)

        # some checks
        solType = t.getType()
        if solType != 'phase':
           logging.warning("Soltab type of "+soltab._v_name+" is: "+solType+" should be phase. Ignoring.")
           continue

        # this will make a selection for the getValues() and getValuesIter()
        userSel = {}
        for axis in t.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        t.setSelection(**userSel)

        names=t.getAxesNames()

        # Collect station properties
        station_dict = H.getAnt(solsetname)
        stations = t.getAxisValues('ant')
        station_positions = np.zeros((len(stations), 3), dtype=np.float)
        for i, station_name in enumerate(stations):
            station_positions[i, 0] = station_dict[station_name][0]
            station_positions[i, 1] = station_dict[station_name][1]
            station_positions[i, 2] = station_dict[station_name][2]
            
        returnAxes=['ant','freq','pol','time']
        for vals, flags, coord, selection in t.getValuesIter(returnAxes=returnAxes,weight=True):

            logging.debug('Flags '+str(np.sum(flags))+' '+str(vals.shape))

            if len(coord['ant']) < 2:
                logging.error('Clock/TEC separation needs at least 2 antennas selected.')
                return 1
            if len(coord['freq']) < 10:
                logging.error('Clock/TEC separation needs at least 10 frequency channels, preferably distributed over a wide range')
                return 1

            freqs=coord['freq']
            stations=coord['ant']
            times=coord['time']

            axes=[i for i in names if i in returnAxes]
            clock,tec,offset,newstations=doFit(vals,np.logical_not(flags),freqs,stations,station_positions,axes,\
                    flagBadChannels=flagBadChannels,flagcut=flagCut,chi2cut=chi2cut,combine_pol=combinePol,doFitoffset=fitOffset,removePhaseWraps=removePhaseWraps)
            weights=tec>-5
            tec[np.logical_not(weights)]=0
            clock[np.logical_not(weights)]=0

            # if pol are combined here we duplicate the output and save same values in XX and YY
            if combinePol:
                #print "clock",clock.shape
                #print "vals",vals.shape
                clock = np.tile(clock,(1,1,2))
                tec = np.tile(tec,(1,1,2))
                offset = np.tile(offset,(1,2))
                weights = np.tile(weights,(1,1,2))

            tf_st = H.makeSoltab(solsetname, 'tec',
                                 axesNames=['time', 'ant','pol'], axesVals=[times, newstations, ['XX','YY']],
                                 vals=tec,
                                 weights=weights)
            sw = solWriter(tf_st)
            sw.addHistory('CREATE (by CLOCKTECFIT operation)')
            tf_st = H.makeSoltab(solsetname, 'clock',
                                 axesNames=['time', 'ant','pol'], axesVals=[times, newstations, ['XX','YY']],
                                 vals=clock*1e-9,
                                 weights=weights)
            sw = solWriter(tf_st)
            sw.addHistory('CREATE (by CLOCKTECFIT operation)')
            tf_st = H.makeSoltab(solsetname, 'phase_offset',
                                 axesNames=['ant','pol'], axesVals=[newstations, ['XX','YY']],
                                 vals=offset,
                                 weights=np.ones_like(offset))
            sw = solWriter(tf_st)
            sw.addHistory('CREATE (by CLOCKTECFIT operation)')

    return 0
