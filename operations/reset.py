#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation reset all the selected amplitudes to 1
# and all other selected solution types to 0

from operations_lib import *
import logging
from h5parm import solWriter

logging.debug('Loading RESET module.')

def run( step, parset, H ):

   soltabs = getParSoltabs( step, parset, H )
   ants = getParAxis( step, parset, H, 'ant' )
   pols = getParAxis( step, parset, H, 'pol' )
   dirs = getParAxis( step, parset, H, 'dir' )

   setWeight = parset.getBool('.'.join(["LoSoTo.Steps", step, "Weight"]), False )

   for soltab in openSoltabs( H, soltabs ):
        t = solWriter(soltab)
        t.setSelection(ant=ants, pol=pols, dir=dirs)
        solType = t.getType()

        if setWeight:
            t.setValues(1., weight=setWeight)
        elif solType == 'amplitude':
            t.setValues(1.)
        else:
            t.setValues(0.)

        t.addHistory('RESET')
   return 0
