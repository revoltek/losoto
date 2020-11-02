#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Fill a solution table with beam response

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading LOFARBEAM module.')

# this funct is called by losoto to set parameters and call the real run()
def _run_parser(soltab, parser, step):
    ms = parser.getstr( step, 'ms') # no default
    inverse = parser.getbool( step, 'inverse', False)
    useElementResponse = parser.getbool( step, 'useElementResponse', True )
    useArrayFactor = parser.getbool( step, 'useArrayFactor', True )
    useChanFreq = parser.getbool( step, 'useChanFreq', True )

    parser.checkSpelling( step, soltab, ['ms', 'inverse', 'useElementResponse', 'useArrayFactor', 'useChanFreq'])
    return run(soltab, ms, inverse, useElementResponse, useArrayFactor, useChanFreq)

def run( soltab, ms, inverse=False, useElementResponse=True, useArrayFactor=True, useChanFreq=True ):
    """
    Generic unspecified step for easy expansion.

    Parameters
    ----------
    opt1 : float
        Is a mandatory parameter.

    """
    # load specific libs
    import numpy as np
    import casacore.tables as pt
    from lofar.stationresponse import stationresponse

    sr = stationresponse(ms, inverse, useElementResponse, useArrayFactor, useChanFreq)

    numants = pt.taql('select gcount(*) as numants from '+ms+'::ANTENNA').getcol('numants')[0]
    times = soltab.getAxisValues('time')

    for vals, coord, selection in soltab.getValuesIter(returnAxes=['ant','time','pol','freq'], weight=False):
        vals = reorderAxes( vals, soltab.getAxesNames(), ['ant','time','freq','pol'])

        for stationnum in range(numants):
            logging.debug('Working on station number %i' % stationnum)
            for itime, time in enumerate(times):
                beam = sr.evaluateStation(time=time, station=stationnum)
                # Reshape from [nfreq, 2, 2] to [nfreq, 4]
                beam = beam.reshape(beam.shape[0], 4)

                if soltab.getAxisLen('pol') == 2:
                    beam = beam[:,[0,3]] # get only XX and YY

                if soltab.getType() == 'amplitude':
                    vals[stationnum, itime, :, :] = np.abs(beam)
                elif soltab.getType() == 'phase':
                    vals[stationnum, itime, :, :] = np.angle(beam)
                else:
                    logging.error('Beam prediction work only for amplitude/phase solution tables.')
                    return 1

        vals = reorderAxes( vals, ['ant','time','freq','pol'], [ax for ax in soltab.getAxesNames() if ax in ['ant','time','freq','pol']] )
        soltab.setValues(vals, selection)

    return 0
 
 
