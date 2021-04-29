#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading DELETEAXIS module.')

def _run_parser(soltab, parser, step):
    axisDelete = parser.getstr( step, 'axisDelete' ) # no default
    fromCell = parser.getstr( step, 'fromCell', '' )

    parser.checkSpelling( step, soltab, ['axisDelete', 'fromCell'])
    return run(soltab, axisDelete, fromCell)

def run( soltab, axisDelete, fromCell):
    """
    Delete an axis only keeping the values of a certain slice.

    Parameters
    ----------
    axisDelete : str
        Axis to delete.

    fromCell : str
        A cell value in axisDelete from which to keep the data values. If it is the string
        "first"/"last" then uses the first/last element of the axis.
    """
    import numpy as np

    if not axisDelete in soltab.getAxesNames():
        logging.error('Cannot find axis %s.' % axisDelete)
        return 1

    if fromCell == 'first':
        fromCell = soltab.getAxisValues(axisDelete)[0]
    elif fromCell == 'last':
        fromCell = soltab.getAxisValues(axisDelete)[-1]

    axisType = type(soltab.getAxisValues(axisDelete)[0])
    try:
        fromCell = np.array([fromCell]).astype(axisType)[0]
    except:
        logging.error('Cannot convert to type %s the value in fromCell: %s.' % (str(axisType),fromCell))
        return 1

    if not fromCell in soltab.getAxisValues(axisDelete):
        logging.error('Cannot find %s in %s.' % (fromCell, axisDelete))
        return 1

    logging.info("Delete axis on soltab: "+soltab.name)

    # get slice with 1 value to replicate
    soltab.setSelection(**{axisDelete:fromCell})
    axisDeleteIdx = soltab.getAxesNames().index(axisDelete)
    fromCellIdx = list(soltab.getAxisValues(axisDelete)).index(fromCell)
    vals = soltab.getValues(retAxesVals=False)
    vals = np.take(soltab.getValues(retAxesVals=False), fromCellIdx, axisDeleteIdx)
    weight = np.take(soltab.getValues(retAxesVals=False, weight=True), fromCellIdx, axisDeleteIdx)
    solset = soltab.getSolset()
    sttype = soltab.getType()
    stname = soltab.name

    axes = soltab.getAxesNames()
    axes.remove(axisDelete)
    axesVals = [soltab.getAxisValues(ax) for ax in axes]
    soltab.delete()

    st = solset.makeSoltab(sttype, stname, axesNames=axes, axesVals=axesVals, vals=vals, weights=weight)

    st.addHistory('DELETEAXIS (over axis %s from cell %s)' % (axisDelete, str(fromCell)))
    return 0
