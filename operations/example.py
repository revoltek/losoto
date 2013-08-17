#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an example operation for LoSoTo

import numpy as np
import logging
from operations_lib import *
from h5parm import solFetcher

logging.info('Loading EXAMPLE module.')

def run( step, parset, H ):
   """
   Generic unspecified step for easy expansion.
   """
   # all the following are LoSoTo function to extract information from the parset

   # get involved solsets using local step values or global values or all
   solsets = getParSolsets( step, parset, H )
   logging.info('Solset: '+str(solsets))
   # get involved soltabs using local step values or global values or all
   soltabs = getParSoltabs( step, parset, H )
   logging.info('Soltab: '+str(soltabs))
   # get list of Antennas using local step values or global values or all
   ants = getParAnts( step, parset, H )
   logging.info('Ant: '+str(ants))
   # get list of Polarizations using local step values or global values or all
   pols = getParPols( step, parset, H )
   logging.info('Pol: '+str(pols))
   # get list of SolTypes using local step values or global values or all
   solTypes = getParSolTypes( step, parset, H )
   logging.info('SolType: '+str(solTypes))
   # get list of Directions using local step values or global values or all
   dirs = getParDirs( step, parset, H )
   logging.info('Dir: '+str(dirs))


   # do something on every soltab (use the openSoltab LoSoTo function)
   for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab.name)
        # use the solFetcher from the H5parm lib
        t = solFetcher(soltab)

        axisNames = t.getAxes()
        logging.info("Axis names are: "+str(axisNames))

        solType = t.getType()
        logging.info("Soltab type is: "+solType)

        # this can be manually done with e.g. setSelection('(ant=='CS001HBA') & (pol=='XX')')
        t.makeSelection(ant=ants, pol=pols, dir=dirs)
        logging.info("Selection is: "+t.selection)

        # find all axis values
        logging.info("Antennas are: "+str(t.getValuesAxis('ant')))
        # but one can also use
        logging.info("Antennas (other method) are: "+str(t.ant))
        logging.info("Frequencies are: "+str(t.freq))
        logging.info("Directions are: "+str(t.dir))
        logging.info("Polarizations are: "+str(t.pol))
        # try to access a non-existent axis
        t.getValuesAxis('nonexistantaxis')

        # now get all values given this selection
        logging.info("Get data using t.val")
        val = t.val
        logging.info("$ val is "+str(val[100,0,0,0,0]))
        #print len(val)
        flag = t.flag
        #print len(flag)
        time = t.time
        thisTime = t.time[100]
        #print len(time)

        # get a rowIterator given a different selection
        logging.info("Get data using getRowsIterator()")
        ants.append('CS002LBA') # change selection
        t.makeSelection(ant=ants, pol=pols, dir=dirs)
        logging.info("Selection is: "+t.selection)
        for i, row in enumerate(t.getRowsIterator()):
            if row['ant'] == 'CS001LBA' and row['time'] == thisTime:
                logging.info("$ val is "+str(row['val']))
            if row['ant'] == 'CS002LBA' and row['time'] == thisTime:
                # update a specific cell value
                row['val'] = '123456'
                row.update()

        # another way to get the data is using the getValuesGrid()
        logging.info("Get data using getValuesGrid()")
        grid, axis = t.getValuesGrid(selection='') # note we reset the selection
        # axis names
        logging.info("Axes: "+str(t.getAxes()))
        # axis shape
        print [len(i) for i in axis]
        # data array shape (same of axis shape)
        print grid.shape
        logging.info("$ val is "+str(grid[100,0,0,1,1]))

        # reset selection
        t.setSelection('')
        logging.info('Reset selection to \'\'')
        logging.info("Antennas (other method) are: "+str(t.ant))
        logging.info("Frequencies are: "+str(t.freq))
        logging.info("Directions are: "+str(t.dir))
        logging.info("Polarizations are: "+str(t.pol))

        # probably the fastest way to dump all the data
        a=[row.fetch_all_fields() for row in t.t.where('(ant == \'CS002LBA\')')]

        return 0 # if everything went fine, otherwise 1


