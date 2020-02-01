#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading BANDPASSTEC module.')

def _run_parser(soltab, parser, step):
    soltabOutTec = parser.getstr( step, 'soltabOutTEC', 'tec000' )
    soltabOutBP = parser.getstr( step, 'soltabOutBP', 'phase000' )
    refAnt = parser.getstr( step, 'refAnt', '')

    parser.checkSpelling( step, soltab, ['soltabOutTEC', 'soltabOutBP', 'refAnt', 'maxResidual'])
    return run(soltab, soltabOutTEC, soltabOutBP, refAnt)


def run( soltab, soltabOutTEC='tec000', soltabOutBP='phase000', refAnt=''):
    """
    Isolate BP from TEC in a calibrator (uGMRT data)

    Parameters
    ----------
    soltabOutTEC : str, optional
        output TEC table name (same solset), by deault "tec000".

    soltabOutBP : str, optional
        output bandpass table name (same solset), by deault "phase000".

    refAnt : str, optional
        Reference antenna, by default get the closest for each antenna.

    """
    import numpy as np

    logging.info("Find BANDPASS+TEC for soltab: "+soltab.name)

    # input check
    solType = soltab.getType()
    if solType != 'phase':
       logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase. Ignoring.")
       return 1

    ants = soltab.getAxisValues('ant')
    if refAnt != '' and refAnt != 'closest' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.error('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[0]
    if refAnt == '': refAnt = ants[0]

    # create new table
    solset = soltab.getSolset()
    soltaboutTEC = solset.makeSoltab(soltype = 'tec', soltabName = soltabOut, axesNames=['ant','time'], \
                      axesVals=[soltab.getAxisValues(axisName) for axisName in ['ant','time']], \
                      vals=np.zeros(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('time'))), \
                      weights=np.ones(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('time'))) )
    soltaboutTEC.addHistory('Created by BANDPASSTEC operation from %s.' % soltab.name)

    soltaboutBP = solset.makeSoltab(soltype = 'phase', soltabName = soltabOut, axesNames=['ant','freq','pol'], \
                      axesVals=[soltab.getAxisValues(axisName) for axisName in ['ant','freq', 'pol']], \
                      vals=np.zeros(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('freq'),soltab.getAxisLen('pol'))), \
                      weights=np.ones(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('freq'),soltab.getAxisLen('pol'))) )
    soltaboutBP.addHistory('Created by BANDPASSTEC operation from %s.' % soltab.name)
 
    # get values
    vals, axesVals = soltab.getValues(retAxesVals=True) 

    # separate bandpass/tec



    plot = False
    if plot:
        pass

    # write solutions back
    soltaboutTEC.setValues( fitd )
    soltaboutBP.setValues( fitd )
    soltaboutTEC.setValues( fitweights, weight=True )
    soltaboutBP.setValues( fitweights, weight=True )

    return 0
