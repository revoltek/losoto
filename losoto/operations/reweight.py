#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation reset the weight vals
# So far it sets weights to a specified number, in the future an algorithm to
# weight solutions properly should be implemented

from losoto.operations_lib import *
import logging

logging.debug('Loading REWEIGHT module.')

def run( step, parset, H ):

   from losoto.h5parm import solWriter
   import numpy as np

   soltabs = getParSoltabs( step, parset, H )

   weightVal = parset.getFloat('.'.join(["LoSoTo.Steps", step, "WeightVal"]), 1. )
   mergeSoltab = parset.getString('.'.join(["LoSoTo.Steps", step, "MergeFromSoltab"]), '' )
   flagBad = parset.getBool('.'.join(["LoSoTo.Steps", step, "FlagBad"]), False )

   for soltab in openSoltabs( H, soltabs ):

        logging.info("Reweighting soltab: "+soltab._v_name)

        sw = solWriter(soltab)

        # axis selection
        userSel = {}
        for axis in sw.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sw.setSelection(**userSel)

        if mergeSoltab != '':
            mss, mst = mergeSoltab.split('/')
            msf = solFetcher(H.getSoltab(mss, mst))
            msf.setSelection(**userSel)
            sf = solFetcher(soltab)
            sf.setSelection(**userSel)
            weights, axes = sf.getValues(weight = True)
            mergeWeights, mergeAxes = msf.getValues(weight = True)
            if axes.keys() != mergeAxes.keys() or weights.shape != mergeWeights.shape:
                logging.error('Impossible merge two tables with different axes values')
                return 1
            weights[ np.where(mergeWeights == 0) ] = 0.
            sw.addHistory('WEIGHT merged from '+mergeSoltab+' for selection:'+str(userSel))
        else:
            weights = weightVal
            sw.addHistory('REWEIGHTED to '+str(weightVal)+' for selection:'+str(userSel))

        sw.setValues(weights, weight=True)

        if flagBad:
            sf = solFetcher(soltab)
            sf.setSelection(**userSel)
            weights = sf.getValues(weight = True, retAxesVals = False)
            vals = sf.getValues(retAxesVals = False)
            if sf.getType() == 'amplitude': weights[np.where(vals == 1)] = 0
            else: weights[np.where(vals == 0)] = 0
            sw.setValues(weights, weight=True)

   return 0
