#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is a script for copying tables values in LoSoTo

import logging
from operations_lib import *

logging.debug('Loading COPY module.')

def run( step, parset, H ):
    """
    Copy values from a table to another (of the same kind)
    If tables have different sampling, then resample the values
    """
    import numpy as np
    from h5parm import solFetcher, solWriter
    
    inTable = parset.getString('.'.join(["LoSoTo.Steps", step, "InTable"]), '' )
    outTable = parset.getString('.'.join(["LoSoTo.Steps", step, "OutTable"]), '' )

    tr = solFetcher(inTable)
    inVals, inCoord = tr.getValues()
    tw = solWriter(outTable)
    outVals, outCoord = tw.getValues()

    shapeDiff = tuple(np.array(outVals.shape)/np.array(inVals.shape))

    inValsNew = np.kron(inVals, np.ones(shapeDiff))
    if inValsNew.shape != outVals.shape:
        logging.error("Incompatible shapes for copying solutions. The outTable axes can be a multiple of the inTable axis dimension.")

    # writing on the output table
    tw.setValues(inValsNew)

    tw.addHistory('COPY (from table %s)' % (inTable))
    return 0
