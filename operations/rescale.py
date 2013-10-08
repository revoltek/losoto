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
    import itertools
    import scipy.interpolate
    import numpy as np
    from h5parm import solFetcher, solWriter

    solsets = getParSolsets( step, parset, H )
    soltabs = getParSoltabs( step, parset, H )
    ants = getParAnts( step, parset, H )
    pols = getParPols( step, parset, H )
    solTypes = getParSolTypes( step, parset, H )
    dirs = getParDirs( step, parset, H )

    calSoltab = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "CalSoltab"]), [] )
    medAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "MedAxes"]), ['time'] )
    interpAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "InterpAxes"]), ['freq'] )
    interpMethod = parset.getString('.'.join(["LoSoTo.Steps", step, "InterpMethod"]), 'linear' )
    if interpMethod not in ["nearest", "linear", "cubic"]:
        logging.error('Interpolation method must be nearest, linear or cubic.')
        return 1

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab.name)

        tr = solFetcher(soltab)
        tw = solWriter(soltab)
        css, cst = calSoltab.split('/')
        cr = solFetcher(H.getSoltab(css, cst))

        axesNames = tr.getAxes()
        for interpAxis in interpAxes:
            if interpAxis not in axesNames:
                logging.error('Axis '+interpAxis+' not found.')
                return 1
        for avgAxis in avgAxes:
            if avgAxis not in axesNames:
                logging.error('Axis '+avgAxis+' not found.')
                return 1


        tr.makeSelection(ant=ants, pol=pols, dir=dirs)
        for vals, coord, nrows in tr.getIterValuesGrid(returnAxes=medAxes+interpAxes, return_nrows=True):

            # constract grid
            coordSel = removeKeys(coord, interpAxes)
            logging.debug("Working on coords:"+str(coordSel))
            cr.makeSelection(**coordSel)
            calValues, calCoord = cr.getValuesGrid()
            #med = np.median(calValues[....])
            print "Med:",med

            # create calibrator coordinates array
            calPoints = []
            for interpAxis in interpAxes:
                calPoints.append(calCoord[interpAxis])
            calPoints = np.array([x for x in itertools.product(*calPoints)])
            # create target coordinates array
            targetPoints = []
            for interpAxis in interpAxes:
                targetPoints.append(coord[interpAxis])
            targetPoints = np.array([x for x in itertools.product(*targetPoints)])

            # interpolation
            valsnew = scipy.interpolate.griddata(calPoints, calValues, targetPoints, interpMethod)
            # fill values outside boudaries with "nearest" solutions
            if interpMethod != 'nearest':
                valsnewNearest = scipy.interpolate.griddata(calPoints, calValues, targetPoints, 'nearest')
                valsnew[ np.where(valsnew == np.nan) ] = valsnewNearest [ np.where(valsnew == np.nan) ]
            valsnew = valsnew.reshape(vals.shape)

            # writing back the solutions
            tw.setValuesGrid(valsnew, nrows)
        tw.flush()

    selection = tw.selection
    tw.addHistory('RESCALE (from table %s with selection %s)' % (calSoltab, selection))
    return 0
