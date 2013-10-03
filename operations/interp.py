#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an interpolation script for LoSoTo

import numpy as np
import logging
import scipy.interpolate
from operations_lib import *
from h5parm import solFetcher

logging.debug('Loading INTERP module.')

def run( step, parset, H ):
    """
    Interpolate the solutions from one table into a destination table
    """
    solsets = getParSolsets( step, parset, H )
    soltabs = getParSoltabs( step, parset, H )
    ants = getParAnts( step, parset, H )
    pols = getParPols( step, parset, H )
    solTypes = getParSolTypes( step, parset, H )
    dirs = getParDirs( step, parset, H )

    calSoltab = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "CalSoltab"]), [] )
    interpAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "InterpAxes"]), ['time','freq'] )
    interpMethod = parset.getString('.'.join(["LoSoTo.Steps", step, "InterpMethod"]), 'linear' )
    if interpMethod not in ["nearest", "linear", "cubic"]:
        logging.error('Interpolation method must be nearest, linear or cubic.')
        return 1

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab.name)

        tr = solFetcher(soltab)
        tw = solWriter(soltab)
        cr = solFetcher(calSoltab)

        axesNames = tr.getAxes()
        for interpAxis in interpAxes:
            if interpAxis not in axesNames:
                logging.error('Axis '+interpAxis+' not found.')
                return 1

        tr.makeSelection(ant=ants, pol=pols, dir=dirs)
        for vals, coord, nrows in tr.getIterValuesGrid(returnAxes=interpAxes, return_nrows=True):

            # constract grid
            cr.makeSelection(**coord)
            calValues, calCoord = cr.getValuesGrid()
            print calValues.shape

            # create calibrator coordinates array
            calPoints = []
            for interpAxis in interpAxes:
                calPoints.append(calCoord[interpAxis])
            # create target coordinates array
            targetPoints = []
            for interpAxis in interpAxes:
                targetPoints.append(coord[interpAxis])

            # interpolation
            valsnew = scipy.interpolate.griddata(calPoints, calValues, targetPoints, interpMethod)
            # fill values outside boudaries with "nearest" solutions
            if interpMethod != 'nearest':
                valsnewNearest = scipy.interpolate.griddata(calPoints, calValues, targetPoints, 'nearest')
                valsnew[ np.where(valsnew == np.nan) ] = valsnewNearest [ np.where(valsnew == np.nan) ]

            # writing back the solutions
            tw.setValuesGrid(valsnew, nrows)
        tw.flush()

    selection = tw.selection
    tw.addHistory('INTERP (from table %s with selection %s)' % (calSoltab, selection))

    return 0
