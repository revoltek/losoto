#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading TECSINGLEFREQ module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'tec000' )
    refAnt = parser.getstr( step, 'refAnt', '')

    parser.checkSpelling( step, soltab, ['soltabOut', 'refAnt'])
    return run(soltab, soltabOut, refAnt)


def run( soltab, soltabOut='tec000', refAnt='' ):
    """
    Estimate TEC from ph solutions assuming no wrap for solution at t=0

    Parameters
    ----------
    soltabOut : str, optional
        output table name (same solset), by deault "tec".

    refAnt : str, optional
        Reference antenna, by default the first.

    """

    import numpy as np
    import scipy.optimize
    from losoto.lib_unwrap import unwrap

    logging.info("Find TEC for soltab: "+soltab.name)

    # input check
    solType = soltab.getType()
    if solType != 'phase':
       logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase. Ignoring.")
       return 1

    ants = soltab.getAxisValues('ant')
    if refAnt != '' and refAnt != 'closest' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.warning('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[0]
    if refAnt == '': refAnt = ants[0]

    # times and ants needs to be complete or selection is much slower
    times = soltab.getAxisValues('time')

    # create new table
    solset = soltab.getSolset()
    soltabout = solset.makeSoltab(soltype = 'tec', soltabName = soltabOut, axesNames=['ant','time','dir'], \
                      axesVals=[soltab.getAxisValues(axisName) for axisName in ['ant','time','dir']], \
                      vals=np.zeros(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('time'),soltab.getAxisLen('dir'))), \
                      weights=np.ones(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('time'),soltab.getAxisLen('dir'))) )
    soltabout.addHistory('Created by TEC operation from %s.' % soltab.name)
        
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['freq','time'], weight=True, refAnt=refAnt):

        assert len(coord['freq']) == 1 # it works with phase at only 1 freq

        if not coord['ant'] == refAnt:

            if (weights == 0.).all() == True:
                logging.warning('Skipping flagged antenna: '+coord['ant'])
                fitweights[:] = 0
            else:

                # unwrap
                vals = np.reshape(unwrap(np.squeeze(vals)), vals.shape)
                vals *= coord['freq']/(-8.44797245e9)
                logging.info('%s: average tec: %f TECU' % (coord['ant'], np.mean(vals)))

        # reorder axes back to the original order, needed for setValues
        soltabout.setSelection(ant=coord['ant'], dir=coord['dir'])
        soltabout.setValues( vals )
        soltabout.setValues( weights, weight=True )

    return 0
