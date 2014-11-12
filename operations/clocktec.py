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
    logging.info('Soltab: '+str(soltabs))
     
    # do something on every soltab (use the openSoltab LoSoTo function)
    #for soltab in openSoltabs( H, soltabs ):
    for soltabname in soltabs:
        solsetname=soltabname.split('/')[0]
        soltab=H.getSoltab(solset=solsetname, soltab=soltabname.split('/')[1])
        logging.info("--> Working on soltab: "+soltab._v_name)
        # use the solFetcher from the H5parm lib
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
        logging.info("Selection is: "+str(t.selection))

        names=t.getAxesNames()
        logging.info("Axis names"+str(t.getAxesNames()))

        # Collect station properties
        station_dict = H.getAnt(solsetname)
        stations = t.getAxisValues('ant')
        station_positions = np.zeros((len(stations), 3), dtype=np.float)
        for i, station_name in enumerate(stations):
            station_positions[i, 0] = station_dict[station_name][0]
            station_positions[i, 1] = station_dict[station_name][1]
            station_positions[i, 2] = station_dict[station_name][2]
            
        returnAxes=['ant','freq','pol','time']
        for vals, coord in t.getValuesIter(returnAxes=returnAxes):
            freqs=coord['freq']
            stations=coord['ant']
            times=coord['time']
            ph=vals[:]
            axes=[i for i in names if i in returnAxes]
            clock,tec,offset,newstations=doFit(ph,freqs,stations,station_positions,axes)
            tf_st = H.makeSoltab(solsetname, 'tec',
                                 axesNames=['time', 'ant','pol'], axesVals=[times, newstations, ['XX','YY']],
                                 vals=tec,
                                 weights=np.ones_like(tec))
            sw = solWriter(tf_st)
            sw.addHistory('CREATE (by CLOCKTECFIT operation)')
            tf_st = H.makeSoltab(solsetname, 'clock',
                                 axesNames=['time', 'ant','pol'], axesVals=[times, newstations, ['XX','YY']],
                                 vals=clock,
                                 weights=np.ones_like(clock))
            sw = solWriter(tf_st)
            sw.addHistory('CREATE (by CLOCKTECFIT operation)')
            tf_st = H.makeSoltab(solsetname, 'phase_offset',
                                 axesNames=['ant','pol'], axesVals=[newstations, ['XX','YY']],
                                 vals=offset,
                                 weights=np.ones_like(offset))
            sw = solWriter(tf_st)
            sw.addHistory('CREATE (by CLOCKTECFIT operation)')

    return 0 # if everything went fine, otherwise 1


