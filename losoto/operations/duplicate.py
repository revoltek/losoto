#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading DUPLICATE module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', '' )
    overwrite = parser.getbool( step, 'overwrite', False )

    parser.checkSpelling( step, soltab, ['soltabOut','overwrite'])
    return run(soltab, soltabOut, overwrite)

def run( soltab, soltabOut='', overwrite=False):
    """
    Duplicate a table

    Parameters
    ----------
    soltabOut : str, optional
        Output table name. By default choose next available from table type.
    overwrite : bool, optional
        Overwrite soltabOut if it already exists?
    """

    if soltabOut == '':
        soltabOut = None

    solset = soltab.getSolset()
    if soltabOut in solset.getSoltabNames() and overwrite:
        logging.info('Overwriting soltabOut {}'.format(soltabOut))
        solset.getSoltab(soltabOut).delete()

    soltabout = solset.makeSoltab(soltype = soltab.getType(), soltabName = soltabOut, axesNames=soltab.getAxesNames(),
        axesVals=[soltab.getAxisValues(axisName) for axisName in soltab.getAxesNames()],
        vals=soltab.getValues(retAxesVals = False), weights=soltab.getValues(weight = True, retAxesVals = False))
    # parmdbType=soltab.obj._v_attrs['parmdb_type'] # deprecated

    logging.info('Duplicate %s -> %s' % (soltab.name, soltabout.name) )

    soltabout.addHistory('DUPLICATE from table %s' % (soltab.name))
    return 0
