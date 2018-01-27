#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an interpolation script for LoSoTo
# WEIGHT: not ready

import logging
from losoto.lib_operations import *

logging.debug('Loading INTERP module.')

def run( step, parset, H ):
    """
    Interpolate the solutions from one table into a destination table
    """
    import itertools
    import scipy.interpolate
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    solTypes = getParSolTypes( step, parset, H )

    calSoltab = parset.getString('.'.join(["LoSoTo.Steps", step, "CalSoltab"]), '' )
    calDir = parset.getString('.'.join(["LoSoTo.Steps", step, "CalDir"]), '' )
    interpAxes = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "InterpAxes"]), ['time','freq'] )
    interpMethod = parset.getString('.'.join(["LoSoTo.Steps", step, "InterpMethod"]), 'linear' )
    medAxis = parset.getString('.'.join(["LoSoTo.Steps", step, "MedAxis"]), '' )
    rescale = parset.getBool('.'.join(["LoSoTo.Steps", step, "Rescale"]), False )

    if interpMethod not in ["nearest", "linear", "cubic"]:
        logging.error('Interpolation method must be nearest, linear or cubic.')
        return 1

    if rescale and medAxis == '':
        logging.error('A medAxis is needed for rescaling.')
        return 1

    # open calibration table
    css, cst = calSoltab.split('/')
    cr = solFetcher(H.getSoltab(css, cst))
    cAxesNames = cr.getAxesNames()

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Interpolating soltab: "+soltab._v_name)

        tr = solFetcher(soltab)
        tw = solWriter(soltab)

        axesNames = tr.getAxesNames()
        for i, interpAxis in enumerate(interpAxes[:]):
            if interpAxis not in axesNames or interpAxis not in cAxesNames:
                logging.error('Axis '+interpAxis+' not found. Ignoring.')
                del interpAxes[i]
        if rescale and (medAxis not in axesNames or medAxis not in cAxesNames):
            logging.error('Axis '+medAxis+' not found. Cannot proceed.')
            return 1

        # axis selection
        userSel = {}
        for axis in tr.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        tr.setSelection(**userSel)

        for vals, coord, selection in tr.getValuesIter(returnAxes=interpAxes):

            # construct grid
            coordSel = removeKeys(coord, interpAxes)
            logging.debug("Working on coords:"+str(coordSel))
            # change dir if sepcified
            if calDir != '':
                coordSel['dir'] = calDir
            cr.setSelection(**coordSel)
            calValues, calCoord = cr.getValues()

            # fill medAxis with the median value
            if rescale:
                axis = cAxesNames.index(medAxis)
                calValues = np.repeat( np.expand_dims( np.median( calValues, axis ), axis ), calValues.shape[axis], axis )

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
            valsNew = scipy.interpolate.griddata(calPoints, calValues, targetPoints, interpMethod)

            # fill values outside boudaries with "nearest" solutions
            # NOTE: in 1D is useless due to scipy bug but this is compensated in the next if
            if interpMethod != 'nearest':
                valsNewNearest = scipy.interpolate.griddata(calPoints, calValues, targetPoints, 'nearest')
                # NaN is != from itself
                valsNew[ np.isnan(valsNew) ] = valsNewNearest [ np.isnan(valsNew) ]

            # fix bug in Scipy which put NaNs outside the convex hull in 1D for 'nearest'
            if len(np.squeeze(vals).shape) == 1:
                import scipy.cluster.vq as vq
                valsNew = np.squeeze(valsNew)
                code, dist = vq.vq(targetPoints, calPoints)
                # NaN is != from itself
                valsNew[ np.isnan(valsNew) ] = calValues[code][ np.isnan(valsNew) ]
                valsNew = valsNew.reshape(vals.shape)

            if rescale:
                # rescale solutions
                axis = interpAxes.index(medAxis)
                valsMed = np.repeat( np.expand_dims( np.median( vals, axis ), axis ), vals.shape[axis], axis )
                valsNewMed = np.repeat( np.expand_dims( np.median( valsNew, axis ), axis ), valsNew.shape[axis], axis )
                valsNew = vals*valsNewMed/valsMed
                #print "Rescaling by: ", valsNewMed[:,0]/valsMed[:,0]

            # writing back the solutions
            tw.selection = selection
            tw.setValues(valsNew)

    tw.addHistory('INTERP (from table %s)' % (calSoltab))
    return 0
