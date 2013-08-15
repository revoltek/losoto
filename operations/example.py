#!/usr/bin/env python
# -*- coding: utf-8 -*-

from operations_lib import *
import logging

logging.info('Loading EXAMPLE module.')

def run( step, parset, H ):
   """
   Generic unspecified step for easy expansion.
   """
   # get involved solsets using local step values or global values or all
   solset = getParSolsets( step, parset, H )
   logging.info('Solset: '+solset)
   # get involved soltabs using local step values or global values or all
   soltab = getParSoltabs( step, parset, H )
   # get list of Antennas using local step values or global values or all
   ant = getParAnts( step, parset, H )
   # get list of Polarizations using local step values or global values or all
   pol = getParPols( step, parset, H )
   # get list of SolTypes using local step values or global values or all
   solType = getParSolTypes( step, parset, H )
   # get list of Directions using local step values or global values or all
   dir = getParDir( step, parset, H )

   return 0 # if everything went fine, otherwise 1

