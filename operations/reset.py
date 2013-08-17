#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Set of operations that LoSoTo can perform
# Each operation is a function

from operations_lib import *
import logging
from h5parm import solFetcher

logging.info('Loading RESET module.')

def run( step, parset, H ):
   """Set specific solutions to 1 (amp) and 0 for all the others
   """
   soltabs = getParSoltabs( step, parset, H )
   ant = getParAnts( step, parset, H )
   pol = getParPols( step, parset, H )
   dir = getParDirs( step, parset, H )
    
   for soltab in openSoltabs( H, soltabs ):
        t = solWriter(soltab)
        t.makeSelection(ant=ant, pol=pol, dir=dir)
        if t.getType() == 'amplitude':
            t.setAxes('val', 1.)
        else:
            t.setAxes('val', 0.)

   return 0
