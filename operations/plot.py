#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implements basic plotting
# 

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import logging
from operations_lib import *
from h5parm import solFetcher

logging.debug('Loading PLOT module.')

def run( step, parset, H ):

   soltabs = getParSoltabs( step, parset, H )
   ants = getParAnts( step, parset, H )
   pols = getParPols( step, parset, H )
   dirs = getParDirs( step, parset, H )

   plotType = parset.getString('.'.join(["LoSoTo.Steps", step, "PlotType"]), '' )
   axesToPlot = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), '' )
   minZ, maxZ = parset.getDoubleVector('.'.join(["LoSoTo.Steps", step, "MinMax"]), [0,0] )
   prefix = parset.getString('.'.join(["LoSoTo.Steps", step, "Prefix"]), '' )

   for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        logging.info("Plotting soltab: "+soltab.name)

        sf.makeSelection(ant=ants, pol=pols, dir=dirs)

        # some checks
        for axis in axesToPlot:
            if axis not in sf.getAxes():
                logging.error('Axis \"'+axis+'\" not found.')
                return 1

        if (len(axesToPlot) != 2 and plotType == '2D') or \
           (len(axesToPlot) != 1 and plotType == '1D'):
            logging.error('Wrong number of axes.')
            return 1

        for vals, coord in sf.getIterValuesGrid(returnAxes=axesToPlot):
            # TODO: implement flag control, using different color?
            
            title = ''
            for axis in coord:
                if axis in axesToPlot: continue
                title += str(coord[axis])+'_'
            title = prefix+title[:-1]

            if plotType == '2D':
                fig = plt.figure()
                ax = plt.subplot(111)
                plt.title(title)
                plt.ylabel(axesToPlot[0])
                plt.xlabel(axesToPlot[1])
                p = ax.imshow(coord[axesToPlot[1]], coord[axesToPlot[0]], vals)
                if not (minZ == 0 and maxZ == 0):
                    plt.zlim(zmin=minZ, zmax=maxZ)
                plt.savefig(title+'.png')
                logging.info("Saving "+title+'.png')

            if plotType == '1D':
                fig = plt.figure()
                ax = plt.subplot(111)
                plt.title(title)
                plt.ylabel(sf.getType())
                if not (minZ == 0 and maxZ == 0):
                    plt.ylim(ymin=minZ, ymax=maxZ)
                plt.xlabel(axesToPlot[0])
                p = ax.plot(coord[axesToPlot[0]], vals)
                plt.savefig(title+'.png')
                logging.info("Saving "+title+'.png')

   return 0


