#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement the unflagging procedure
# WEIGHT: flag compliant, no need for weight

import logging
from operations_lib import *
import numpy as np

logging.debug('Loading UNFLAG module.')

def run( step, parset, H ):

    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    unflagzeros = parset.getBool('.'.join(["LoSoTo.Steps", step, "UnflagZeros"]), False )
    
    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Unflagging soltab: "+soltab._v_name)

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)
        sw.setSelection(**userSel)

        weights = sf.getValues(retAxesVals = False, weight = True)
        weights = np.ones_like(weights)

        if not unflagzeros:
            vals = sf.getValues(retAxesVals = False)
            if sf.getType() == 'amplitude': weights[np.where(vals == 1)] = 0
            else: weights[np.where(vals == 0)] = 0

        logging.debug('Percentage of data flagged: %.3f %%' % (100.*len(np.where(weights==0)[0])/len(weights.flat)))

        sw.setValues(weights, weight=True)

        sw.addHistory('UNFLAG')
    return 0
