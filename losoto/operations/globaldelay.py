#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading GLOBALDELAY module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'delay000' )
    refAnt = parser.getstr( step, 'refAnt', '')

    parser.checkSpelling( step, soltab, ['soltabOut', 'refAnt'])
    return run(soltab, soltabOut, refAnt)


def run( soltab, soltabOut='tec000', refAnt='' ):
    """
    Bruteforce TEC extraction from phase solutions.

    Parameters
    ----------
    soltabOut : str, optional
        output table name (same solset), by deault "tec".

    refAnt : str, optional
        Reference antenna, by default the first.

    """
    import numpy as np
    import scipy.optimize

    def mod(d):
        return np.mod(d + np.pi, 2.*np.pi) - np.pi

    def cost_f(d, freq, y):
        nfreq, ntime = y.shape
        phase = mod(2*np.pi*d*freq).repeat(ntime).reshape(nfreq,ntime)
        dist = np.abs(mod(phase - y))
        ngood = np.sum(~np.isnan(dist))
        return np.nansum(dist/ngood)

    logging.info("Find global DELAY for soltab: "+soltab.name)

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
    soltabout = solset.makeSoltab(soltype = 'clock', soltabName = soltabOut, axesNames=['ant','time'], \
                      axesVals=[soltab.getAxisValues(axisName) for axisName in ['ant','time']], \
                      vals=np.zeros(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('time'))), \
                      weights=np.ones(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('time'))) )
    soltabout.addHistory('Created by GLOBALDELAY operation from %s.' % soltab.name)
        
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['freq','time'], weight=True, refAnt=refAnt):

        if len(coord['freq']) < 5:
            logging.error('Delay estimation needs at least 5 frequency channels, preferably distributed over a wide range.')
            return 1

        # reorder axes
        vals = reorderAxes( vals, soltab.getAxesNames(), ['freq','time'] )
        weights = reorderAxes( weights, soltab.getAxesNames(), ['freq','time'] )

        ranges = (-1e-7,1e-7)
        Ns = 1001

        delay_fitresult = np.zeros(len(times))
        weights_fitresult = np.ones(len(times))

        if not coord['ant'] == refAnt:

            if (weights == 0.).all() == True:
                logging.warning('Skipping flagged antenna: '+coord['ant'])
                weights_fitresult[:] = 0
            else:

                freq = np.copy(coord['freq'])

                # brute force
                fit = scipy.optimize.brute(cost_f, ranges=(ranges,), Ns=Ns, args=(freq, vals))
                delay_fitresult[:] = fit

                logging.info('%s: average delay: %f ns' % (coord['ant'], fit*1e9))

        # reorder axes back to the original order, needed for setValues
        soltabout.setSelection(ant=coord['ant'])
        soltabout.setValues( delay_fitresult )
        soltabout.setValues( weights_fitresult, weight=True )

    return 0
