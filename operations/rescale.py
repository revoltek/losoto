#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is a rescaling script for LoSoTo

import numpy as np
import logging
import scipy.interpolate
from operations_lib import *
from h5parm import solFetcher

logging.debug('Loading RESCALE module.')

def run( step, parset, H ):
    """
    Rescale solutions of one table to match the median of another.
    The median of the MedAxes is used to find the rescaling factor.
    The InterpAxes are instead just interpolated.
    """
    solsets = getParSolsets( step, parset, H )
    soltabs = getParSoltabs( step, parset, H )
    ants = getParAnts( step, parset, H )
    pols = getParPols( step, parset, H )
    solTypes = getParSolTypes( step, parset, H )
    dirs = getParDirs( step, parset, H )

    calSoltab = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "CalSoltab"]), [] )
    medAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "MedAxes"]), ['time'] )
    interpAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "InterpAxes"]), ['freq'] )

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab.name)

        t = solFetcher(soltab)
        ct = solFetcher(calSoltab)

        axesNames = t.getAxes()
        for avgAxis in avgAxes:
            if avgAxis not in axesNames:
                logging.error('Axis '+interpAxis+' not found.')
                return 1

        t.makeSelection(ant=ants, pol=pols, dir=dirs)
        for vals, coord, nrows in sf.getIterValuesGrid(returnAxes=avgAxes, return_nrows=True):

            # constract grid
            ct.makeSelection(**coord)
            calValues, calCoord = ct.getValuesGrid()
            print calValues.shape
            print calCoord
            med = np.median(calValues)
            print "Med:",med

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

            # writing back the solutions
            sw.setValuesGrid(valsnew, nrows)

        # a
    return 0
