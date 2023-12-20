#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Fill a solution table with beam response

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading LOFARBEAM module.')

# this funct is called by losoto to set parameters and call the real run()
def _run_parser(soltab, parser, step):
    ms = parser.getstr( step, 'ms') # no default
    useElementResponse = parser.getbool( step, 'useElementResponse', True )
    useArrayFactor = parser.getbool( step, 'useArrayFactor', True )
    allChan = parser.getbool( step, 'allChan', False)

    parser.checkSpelling( step, soltab, ['ms', 'useElementResponse', 'useArrayFactor', 'allChan'])
    return run(soltab, ms, useElementResponse, useArrayFactor, allChan)

def run( soltab, ms, useElementResponse=True, useArrayFactor=True, allChan=False ):
    """
    Generic unspecified step for easy expansion.

    Parameters
    ----------
    ms : string
        A measurement set

    allChan : bool, optional
        Assume that the frequency axes of the MS and the soltab are the same, faster. By default: False.

    """
    # load specific libs
    import numpy as np
    import everybeam as eb

    tel = eb.load_telescope(ms, element_response_model='hamaker')#, inverse, useElementResponse, useArrayFactor, useChanFreq)

    #import casacore.tables as pt
    #numants = pt.taql('select gcount(*) as numants from '+ms+'::ANTENNA').getcol('numants')[0]
    numants = soltab.getAxisLen('ant')
    times = soltab.getAxisValues('time')
    freqs = soltab.getAxisValues('freq')

    if soltab.getType() != 'amplitude' and soltab.getType() != 'phase':
        logging.error('Beam prediction work only for amplitude/phase solution tables.')
        return 1

    for vals, coord, selection in soltab.getValuesIter(returnAxes=['ant','time','pol','freq'], weight=False):
        vals = reorderAxes( vals, soltab.getAxesNames(), ['ant','time','freq','pol'])

        for stationnum in range(numants):
            logging.debug('Working on station number %i' % stationnum)
            for itime, time in enumerate(times):
                if allChan:
                    beam = tel.station_response(time=time, station_idx=stationnum)
                    beam = beam.reshape(beam.shape[0], 4)
                    if soltab.getAxisLen('pol') == 2:
                        beam = beam[:,[0,3]] # get only XX and YY

                    if soltab.getType() == 'amplitude':
                        vals[stationnum, itime, :, :] = np.abs(beam)
                    elif soltab.getType() == 'phase':
                        vals[stationnum, itime, :, :] = np.angle(beam)
                else:
                    for ifreq, freq in enumerate(freqs):
                        beam = tel.station_response(time=time, station_idx=stationnum, freq=freq)
                        # Reshape from [2, 2] to [4]
                        beam = beam.reshape(4)
    
                        if soltab.getAxisLen('pol') == 2:
                            beam = beam[[0,3]] # get only XX and YY
    
                        if soltab.getType() == 'amplitude':
                            vals[stationnum, itime, ifreq, :] = np.abs(beam)
                        elif soltab.getType() == 'phase':
                            vals[stationnum, itime, ifreq, :] = np.angle(beam)


        vals = reorderAxes( vals, ['ant','time','freq','pol'], [ax for ax in soltab.getAxesNames() if ax in ['ant','time','freq','pol']] )
        soltab.setValues(vals, selection)

    return 0
 
 
