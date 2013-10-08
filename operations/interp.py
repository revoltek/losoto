#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an interpolation script for LoSoTo

import logging
from operations_lib import *

logging.debug('Loading INTERP module.')

def run( step, parset, H ):
    """
    Interpolate the solutions from one table into a destination table
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

    calSoltab = parset.getString('.'.join(["LoSoTo.Steps", step, "CalSoltab"]), '' )
    interpAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "InterpAxes"]), ['time','freq'] )
    interpMethod = parset.getString('.'.join(["LoSoTo.Steps", step, "InterpMethod"]), 'linear' )
    medAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "MedAxes"]), [''] )
    rescale = parset.getBool('.'.join(["LoSoTo.Steps", step, "Rescale"]), False )

    if interpMethod not in ["nearest", "linear", "cubic"]:
        logging.error('Interpolation method must be nearest, linear or cubic.')
        return 1
    
#    for i, medAxis in enumerate(medAxes[:]):
#        if medAxis in interpAxes:
#            logging.warning('Axis '+medAxis+' is an interpolation axis, removing it from medAxes list')
#            del medAxes[i]

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab.name)

        tr = solFetcher(soltab)
        tw = solWriter(soltab)
        css, cst = calSoltab.split('/')
        cr = solFetcher(H.getSoltab(css, cst))

        axesNames = tr.getAxes()
        for i, interpAxis in enumerate(interpAxes[:]):
            if interpAxis not in axesNames:
                logging.error('Axis '+interpAxis+' not found. Ignoring.')
                del interpAxes[i]
        for i, medAxis in enumerate(medAxes[:]):
            if medAxis not in axesNames:
                logging.error('Axis '+medAxis+' not found. Ignoring.')
                del medAxes[i]

        tr.makeSelection(ant=ants, pol=pols, dir=dirs)
        for vals, coord, nrows in tr.getIterValuesGrid(returnAxes=interpAxes, return_nrows=True):

            # constract grid
            coordSel = removeKeys(coord, interpAxes)
            logging.debug("Working on coords:"+str(coordSel))
            cr.makeSelection(**coordSel)
            calValues, calCoord = cr.getValuesGrid()
            print calValues
            for medAxis in medAxes:
                axis = cr.getAxes().index(medAxis)
                print np.repeat(np.expand_dims(np.median( calValues, axis), axis), len(cr.getAxisVal(medAxis), axis))
            print calValues
              
            # create a list of values whose coords are calPoints
            calValues = np.ndarray.flatten(calValues)

            # create calibrator/target coordinates arrays
            calPoints = []
            targetPoints = []
            for interpAxis in interpAxes:
                calPoints.append(calCoord[interpAxis])
                targetPoints.append(coord[interpAxis])
            calPoints = np.array([x for x in itertools.product(*calPoints)])
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
    tw.addHistory('INTERP (from table %s with selection %s)' % (calSoltab, selection))
    return 0
