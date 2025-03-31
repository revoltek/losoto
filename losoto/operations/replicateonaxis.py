#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading REPLICATEONAXIS module.')

def _run_parser(soltab, parser, step):
    axisReplicate = parser.getstr( step, 'axisReplicate' ) # no default
    fromCell = parser.getstr( step, 'fromCell', '' )
    updateWeights = parser.getbool( step, 'updateWeights', False )

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
        A cell value in axisReplicate from which to copy the data values. If it is the string
        "first"/"last" then uses the first/last element of the axis. "nonflaggedCS" can be used
        in the case all CS are the same and one wants to expand those to other stations, this solves
        the issue that might happen if the first/last station has NaNs.

    updateWeights : bool
        If False then weights are untoched, if True they are replicated like data (usable only with fromCell: first or last). Default: False.
    """
    import numpy as np

    if not axisReplicate in soltab.getAxesNames():
        logging.error('Cannot find axis %s.' % axisReplicate)
        return 1

    if fromCell == 'first':
        fromCell = soltab.getAxisValues(axisReplicate)[0]
    elif fromCell == 'last':
        fromCell = soltab.getAxisValues(axisReplicate)[-1]
    elif fromCell == 'nonflaggedCS':
        fromCell = 'CS*'
        if updateWeights:
            logging.error('updateWeights must be false with fromCell=nonflaggedCS.')
            return 1

    if fromCell == 'first' or fromCell == 'last':  
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
    axisReplicateLen = soltab.getAxisLen(axisReplicate, ignoreSelection=False) # keep selection into account
    old_selection = np.copy(soltab.selection)

    # get slice with 1 value to replicate
    soltab.setSelection(**{axisReplicate:fromCell})
    vals = soltab.getValues(retAxesVals=False)
    if updateWeights:
        weights = soltab.getValues(retAxesVals=False, weight=True)

    #cellPos = list(soltab.getAxisValues(axisReplicate)).index(fromCell)
    axisReplicatePos = soltab.getAxesNames().index(axisReplicate)
    if fromCell == 'CS*':
        vals = np.nanmean(vals, axis=axisReplicatePos)

    # expand on the right axis
    vals = np.repeat(vals, repeats=axisReplicateLen, axis=axisReplicatePos)

    # write back
    soltab.selection = old_selection
    soltab.setValues(vals)
    if updateWeights:
        soltab.setValues(weights, weight=True)

    soltab.addHistory('REPLICATEONAXIS (over axis %s from cell %s)' % (axisReplicate, str(fromCell)))
    return 0
