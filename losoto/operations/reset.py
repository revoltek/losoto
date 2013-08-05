#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Set of operations that LoSoTo can perform
# Each operation is a function

from ..operations_lib import *

print "Loaded reset!"

def run( step, parset, H ):
   """Set specific solutions to 1 (amp) and 0 for all the others
   """
   solsets = getParSolsets( step, parset, H )
   print solsets
   soltabs = getParSoltabs( step, parset, H )
   print soltabs
   ant = getParAnts( step, parset, H )
   pol = getParPols( step, parset, H )
   dir = getParDirs( step, parset, H )
    
   for soltab in openSoltabs( H, solsets, soltabs):
        t = solFetcher(soltab)
        print formatSelection(ant=ant, pol=pol, dir=dir)
        t.setSelection(formatSelection(ant=ant, pol=pol, dir=dir))
        r = t.getRowsIterator()


   return 0
