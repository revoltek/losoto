#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an example operation for LoSoTo

# For an example on how to write a parallel version of an operation, look at the flag operation

import logging
from losoto.operations_lib import *

logging.debug('Loading EXAMPLE module.')

def run( step, parset, H ):
   """
   Generic unspecified step for easy expansion.
   """
   import numpy as np
   from losoto.h5parm import solFetcher, solWriter
   # all the following are LoSoTo function to extract information from the parset

   # get involved solsets using local step values or global values or all
   solsets = getParSolsets( step, parset, H )
   logging.info('Solset: '+str(solsets))
   # get involved soltabs using local step values or global values or all
   soltabs = getParSoltabs( step, parset, H )
   logging.info('Soltab: '+str(soltabs))
   # get list of SolTypes using local step values or global values or all
   solTypes = getParSolTypes( step, parset, H )
   logging.info('SolType: '+str(solTypes))


   # do something on every soltab (use the openSoltab LoSoTo function)
   for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab._v_name)
        # use the solFetcher from the H5parm lib
        t = solFetcher(soltab)
        tw = solWriter(soltab)

        axisNames = t.getAxesNames()
        logging.info("Axis names are: "+str(axisNames))

        solType = t.getType()
        logging.info("Soltab type is: "+solType)

        # this will make a selection for the getValues() and getValuesIter()
        # interpret every entry in the parset which has an axis name as a selector
        userSel = {}
        for axis in t.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        t.setSelection(**userSel)

        t.setSelection(ant=ants, pol=pols, dir=dirs)
        logging.info("Selection is: "+str(t.selection))

        # find axis values
        logging.info("Antennas (no selection) are: "+str(t.getAxisValues('ant', ignoreSelection=True)))
        logging.info("Antennas (with selection) are: "+str(t.getAxisValues('ant')))
        # but one can also use (selection is active here!)
        logging.info("Antennas (other method) are: "+str(t.ant))
        logging.info("Frequencies are: "+str(t.freq))
        logging.info("Directions are: "+str(t.dir))
        logging.info("Polarizations are: "+str(t.pol))
        # try to access a non-existent axis
        t.getAxisValues('nonexistantaxis')

        # now get all values given this selection
        logging.info("Get data using t.val")
        val = t.val
        logging.debug('shape of val: '+str(t.val.shape))
        logging.info("$ val is "+str(val[0,0,0,0,100]))
        weight = t.weight
        time = t.time
        thisTime = t.time[100]

        # another way to get the data is using the getValues()
        logging.info("Get data using getValues()")
        grid, axes = t.getValues()
        # axis names
        logging.info("Axes: "+str(t.getAxesNames()))
        # axis shape
        print axes
        print [t.getAxisLen(axis) for axis in axes] # not ordered, is a dict!
        # data array shape (same of axis shape)
        logging.info("Shape of values: "+str(grid.shape))
        #logging.info("$ val is "+str(grid[0,0,0,0,100]))

        # reset selection
        t.setSelection()
        logging.info('Reset selection to \'\'')
        logging.info("Antennas are: "+str(t.ant))
        logging.info("Frequencies are: "+str(t.freq))
        logging.info("Directions are: "+str(t.dir))
        logging.info("Polarizations are: "+str(t.pol))

        # finally the getValuesIter allaws to iterate across all possible combinations of a set of axes
        logging.info('Iteration on time/freq')
        for vals, coord, selection in t.getValuesIter(returnAxes=['time','freq']):
            # writing back the solutions
            tw.selection = selection
            tw.setValues(vals)
        logging.info('Iteration on time')
        for vals, coord, selection in t.getValuesIter(returnAxes=['time']):
            # writing back the solutions
            tw.selection = selection
            tw.setValues(vals)   
        logging.info('Iteration on dir after selection to 1 dir')
        t.setSelection(dir='pointing') 
        for vals, coord, selection in t.getValuesIter(returnAxes=['dir']):
            # writing back the solutions
            tw.selection = selection
            tw.setValues(vals)
 
 
    
   return 0 # if everything went fine, otherwise 1


