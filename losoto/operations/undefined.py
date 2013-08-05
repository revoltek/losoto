#!/usr/bin/env python
# -*- coding: utf-8 -*-

def run( step, parset, H ):
   """
   Generic unspecified step for easy expansion.
   """
   solset = getParSolsets( step, parset, H )
   soltab = getParSoltabs( step, parset, H )
   ant = getParAnts( step, parset, H )
   pol = getParPols( step, parset, H )

   raise Exception('Not yet implemented.')

   return 0 # if everything went fine, otherwise 1

