#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.operations_lib import *

logging.debug('Loading DUPLICATE module.')

def run_parser(soltab, parser, step):
    outTab = parser.getstr( step, 'outtab', '' )
    return run(soltab, outTab)


def run( soltab, outTab=''):
    """
    Duplicate a table

    Parameters
    ----------
    outTab : str, optional
        Output table name. By default choose next available from table type.
    """

    if outTab == '':
        outTab = None

    solset = soltab.getSolset()
    soltabout = solset.makeSoltab(soltype = soltab.getType(), soltabName = outTab, axesNames=soltab.getAxesNames(), \
        axesVals=[soltab.getAxisValues(axisName) for axisName in soltab.getAxesNames()], \
        vals=soltab.getValues(retAxesVals = False), weights=soltab.getValues(weight = True, retAxesVals = False), parmdbType=soltab.obj._v_attrs['parmdb_type'])

    logging.info('Duplicate %s -> %s' % (soltab.name, soltabout.name) )

    soltabout.addHistory('DUPLICATE (from table %s)' % (soltab.name))
    return 0
