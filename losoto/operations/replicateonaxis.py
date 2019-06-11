#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

logging.debug('Loading REPLICATEONAXIS module.')

def _run_parser(soltab, parser, step):
    axisReplicate = parser.getstr( step, 'axisReplicate' ) # no default
    fromCell = parser.getstr( step, 'fromCell', '' )
    updateWeights = parser.getbool( step, 'updateWeights', True )

    parser.checkSpelling( step, soltab, ['axisReplicate', 'fromCell', 'updateWeights'])
    return run(soltab, axisReplicate, fromCell, updateWeights)

def run( soltab, axisReplicate, fromCell, updateWeights=True):
    """
    Replace the values along a certain axis taking them from one specific axic cell

    Parameters
    ----------
    axisReplicate : str
        Axis along which replicate the values.

    fromCell : str
        A cell value in axisReplicate from which to copy the data values.

    updateWeights : bool
        If False then weights are untoched, if True they are replicated like data. Default: True.
    """
    import numpy as np

    if not axisReplicate in soltab.getAxesNames():
        logging.error('Cannot find axis %s.' % axisReplicate)
        return 1

    axisType = type(soltab.getAxisValues(axisReplicate)[0])
    try:
        fromCell = np.array([fromCell]).astype(axisType)[0]
    except:
        logging.error('Cannot convert to type %s the value in fromCell: %s.' % (str(axisType),fromCell))
        return 1

    if not fromCell in soltab.getAxisValues(axisReplicate):
        logging.error('Cannot find %s in %s.' % (fromCell, axisReplicate))
        return 1

    logging.info("Replicate axis on soltab: "+soltab.name)

    # get the cell to replicate
    soltab.setSelection(**{axisReplicate:fromCell})
    vals = soltab.getValues(retAxesVals=False)
    if updateWeights:
        weights = soltab.getValues(retAxesVals=False, weight=True)

    cellPos = list(soltab.getAxisValues(axisReplicate)).index(fromCell)
    axisReplicateLen = soltab.getAxisLen(axisReplicate, ignoreSelection=True)
    axisReplicatePos = soltab.getAxesNames().index(axisReplicate)

    # expand on the right axis
    vals = np.repeat(vals, repeats=axisReplicateLen, axis=axisReplicatePos)

    # write back
    soltab.setSelection()
    soltab.setValues(vals)
    if updateWeights:
        soltab.setValues(weights, weight=True)

    soltab.addHistory('REPLICATEONAXIS (over axis %s)' % (axisReplicate))
    return 0
