#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation reset the weight vals
# So far it sets weights to a specified number, in the future an algorithm to
# weight solutions properly should be implemented

from operations_lib import *
import logging

logging.debug('Loading REWEIGHT module.')

def run( step, parset, H ):

   from h5parm import solWriter
   import numpy as np

   soltabs = getParSoltabs( step, parset, H )

   weightVal = parset.getFloat('.'.join(["LoSoTo.Steps", step, "WeightVal"]), 1. )
   mergeSoltab = parset.getString('.'.join(["LoSoTo.Steps", step, "MergeFromSoltab"]), '' )

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
            newWeights, axes = sf.getValues(weight = True)
            mergeWeights, mergeAxes = msf.getValues(weight = True)
            if axes.keys() != mergeAxes.keys() or newWeights.shape != mergeWeights.shape:
                logging.error('Impossible merge two tables with different axes values')
                return 1
            newWeights[ np.where(mergeWeights == 0) ] = 0.
            sw.setValues(newWeights, weight = True)
            sw.addHistory('WEIGHT merged from '+mergeSoltab+' for selection:'+str(userSel))
        else:
            sw.setValues(weightVal, weight=True)
            sw.addHistory('REWEIGHTED to '+str(weightVal)+' for selection:'+str(userSel))
   return 0
