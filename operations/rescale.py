#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is a rescaling script for LoSoTo

import logging
from operations_lib import *

logging.debug('Loading RESCALE module.')

def run( step, parset, H ):
    """
    Rescale solutions of one table to match the median of another.
    The median of the MedAxes is used to find the rescaling factor.
    The InterpAxes are instead just interpolated.
    """
    import numpy as np
    import scipy.interpolate
    from h5parm import solFetcher

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

        tr = solFetcher(soltab)
        tw = solWriter(soltab)
        cr = solFetcher(calSoltab)

        axesNames = t.getAxes()
        for avgAxis in avgAxes:
            if avgAxis not in axesNames:
                logging.error('Axis '+interpAxis+' not found.')
                return 1

        tr.makeSelection(ant=ants, pol=pols, dir=dirs)
        for vals, coord, nrows in tr.getIterValuesGrid(returnAxes=medAxes+interpAxes, return_nrows=True):

            # constract grid
            cr.makeSelection(**coord)
            calValues, calCoord = cr.getValuesGrid()
            print calValues.shape
            print "coords:", calCoord, coord
            #med = np.median(calValues[....])
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
            tw.setValuesGrid(valsnew, nrows)
        tw.flush()

    selection = tw.selection
    tw.addHistory('RESCALE (from table %s with selection %s)' % (calSoltab, selection))


        # a
    return 0
