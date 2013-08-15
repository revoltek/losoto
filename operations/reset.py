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
   solsets = getParSolsets( step, parset, H )
   print "Solsets", solsets
   soltabs = getParSoltabs( step, parset, H )
   print "Soltabs", soltabs
   ant = getParAnts( step, parset, H )
   print ant
   pol = getParPols( step, parset, H )
   print pol
   dir = getParDirs( step, parset, H )
   print dir
    
   for soltab in openSoltabs( H, soltabs ):
        t = solFetcher(soltab)
        print formatSelection(ant=ant, pol=pol, dir=dir)
        t.setSelection(formatSelection(ant=ant, pol=pol, dir=dir))
        r = t.getRowsIterator()


   return 0
