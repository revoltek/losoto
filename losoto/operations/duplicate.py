#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is a script for duplicate a table in LoSoTo

import logging
from losoto.operations_lib import *

logging.debug('Loading DUPLICATE module.')

def run( step, parset, H ):
    """
    Copy values from a table to another (of the same kind)
    If tables have different sampling, then resample the values
    """
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter
    
    inTable = parset.getString('.'.join(["LoSoTo.Steps", step, "InTable"]), '' ) # complete solset/soltab
    outTable = parset.getString('.'.join(["LoSoTo.Steps", step, "OutTable"]), '' ) # complete solset/soltab or ''

    if inTable == '':
        logging.error('InTable is undefined.')
        return 1

    if outTable == '':
        outSolsetName = inTable.split('/')[0]
        outTableName = None
    else:
        outSolsetName = outTable.split('/')[0]
        outTableName = outTable.split('/')[1]

    ss, st = inTable.split('/')
    sf = solFetcher(H.getSoltab(ss, st))

    t = H.makeSoltab(solset = outSolsetName, soltype = sf.getType(), soltab = outTableName, axesNames=sf.getAxesNames(), \
        axesVals=[sf.getAxisValues(axisName) for axisName in sf.getAxesNames()], \
        vals=sf.getValues(retAxesVals = False), weights=sf.getValues(weight = True, retAxesVals = False), parmdbType=sf.t._v_attrs['parmdb_type'])

    sw = solWriter(t)
    sw.addHistory('DUPLICATE (from table %s)' % (inTable))
    return 0
