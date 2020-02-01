#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading SPLITLEAK module.')

def _run_parser(soltab, parser, step):
    soltabOutG = parser.getstr( step, 'soltabOutG' ) # no default
    soltabOutD = parser.getstr( step, 'soltabOutD' ) # no default

    parser.checkSpelling( step, soltab, ['soltabOutG','soltabOutD'])
    return run(soltab, soltabOutG, soltabOutD)

def run( soltab, soltabOutG=None, soltabOutD=None):
    """
    Duplicate a table

    Parameters
    ----------
    soltabOutG : str, optional
        Output table name (diagonal component). By default choose next available from table type.

    soltabOutD : str, optional
        Output table name (leakage component). By default choose next available from table type.
    """

    logging.info('Split leakage tables %s -> %s + %s' % (soltab.name, soltabOutG, soltabOutD) )

    if soltab.getType() != 'amplitude' and soltab.getType() != 'phase':
        logging.error('SPLITLEAK can work only on amplitude/phase soltabs. Found: %s.' % soltab.getType())
        return 1
    if not np.all(soltab.getAxisValues('pol') == ['XX', 'XY', 'YX', 'YY']):
        logging.error('Pol in unusual order or not linear: not implemented.')
        return 1

    solset = soltab.getSolset()

    ### G component
    soltabOutG = solset.makeSoltab(soltype = soltab.getType(), soltabName = soltabOutG, axesNames=soltab.getAxesNames(), \
        axesVals=[soltab.getAxisValues(axisName) for axisName in soltab.getAxesNames()], \
        vals=soltab.getValues(retAxesVals = False), weights=soltab.getValues(weight = True, retAxesVals = False))

    # set offdiag to 0
    soltabOutG.setSelection( pol=['XY','YX'] )
    soltabOutG.setValues( 0. )

    ### D component
    soltabOutD = solset.makeSoltab(soltype = soltab.getType(), soltabName = soltabOutD, axesNames=soltab.getAxesNames(), \
        axesVals=[soltab.getAxisValues(axisName) for axisName in soltab.getAxesNames()], \
        vals=soltab.getValues(retAxesVals = False), weights=soltab.getValues(weight = True, retAxesVals = False))

    # divide offdiag by diag, then set diag to 1 (see Hamaker+ 96, appendix D)
    soltabOutD.setSelection(pol=['XX','YY'])
    valsDiag = np.copy(soltabOutD.getValues( retAxesVals = False ))
    if soltab.getType() == 'amplitude':
        soltabOutD.setValues( 1. )
    if soltab.getType() == 'phase':
        soltabOutD.setValues( 0. )

    soltabOutD.setSelection(pol=['XY','YX'])
    valsOffdiag = soltabOutD.getValues( retAxesVals = False )
    if soltab.getType() == 'amplitude':
        soltabOutD.setValues( valsOffdiag/valsDiag )
    if soltab.getType() == 'phase':
        soltabOutD.setValues( valsOffdiag - valsDiag )

    soltabOutG.addHistory('SPLITLEAK: G component of %s' % (soltab.name))
    soltabOutD.addHistory('SPLITLEAK: D component of %s' % (soltab.name))
    return 0
