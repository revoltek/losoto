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
   ant = getParAnts( step, parset, H )
   pol = getParPols( step, parset, H )
   dir = getParDirs( step, parset, H )

   for soltab in openSoltabs( H, soltabs ):
        t = solWriter(soltab)
        t.makeSelection(ant=ant, pol=pol, dir=dir)
        solType = t.getType()

        if solType == 'amplitude':
            t.setAxis('val', 1.)
        else:
            t.setAxis('val', 0.)

        t.addHistory('RESET')
   return 0
