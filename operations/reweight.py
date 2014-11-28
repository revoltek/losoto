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

   soltabs = getParSoltabs( step, parset, H )

   weightVal = parset.getFloat('.'.join(["LoSoTo.Steps", step, "WeightVal"]), 1. )

   for soltab in openSoltabs( H, soltabs ):

        logging.info("Reweighting soltab: "+soltab._v_name+" (weight = "+str(weightVal)+")")

        sw = solWriter(soltab)

        # axis selection
        userSel = {}
        for axis in sw.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sw.setSelection(**userSel)

        solType = sw.getType()

        sw.setValues(weightVal, weight=True)

        sw.addHistory('REWEIGHTED to '+str(weightVal)+' for selection:'+str(userSel))
   return 0
