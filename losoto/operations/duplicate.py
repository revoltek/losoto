#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

logging.debug('Loading DUPLICATE module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', '' )
    return run(soltab, soltabOut)

def run( soltab, soltabOut=''):
    """
    Duplicate a table

    Parameters
    ----------
    soltabOut : str, optional
        Output table name. By default choose next available from table type.
    """

    if soltabOut == '':
        soltabOut = None

    solset = soltab.getSolset()
    soltabout = solset.makeSoltab(soltype = soltab.getType(), soltabName = soltabOut, axesNames=soltab.getAxesNames(), \
        axesVals=[soltab.getAxisValues(axisName) for axisName in soltab.getAxesNames()], \
        vals=soltab.getValues(retAxesVals = False), weights=soltab.getValues(weight = True, retAxesVals = False))
    # parmdbType=soltab.obj._v_attrs['parmdb_type'] # deprecated

    logging.info('Duplicate %s -> %s' % (soltab.name, soltabout.name) )

    soltabout.addHistory('DUPLICATE (from table %s)' % (soltab.name))
    return 0
