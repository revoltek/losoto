#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.operations_lib import *
import logging

logging.debug('Loading REWEIGHT module.')

def run_parser(soltab, parser, step):
    weightVal = parser.getfloat( step, 'weightVal', 1. )
    soltabImport = parser.getstr( step, 'soltabImport', '' )
    flagBad = parser.getbool( step, 'flagBad', False )
    return run(soltab, weightVal, soltabImport, flagBad)


def run( soltab, weightVal=1., soltabImport='', flagBad=False ):
    """
    This operation reset the weight vals
    So far it sets weights to a specified number, in the future an algorithm to
    weight solutions properly should be implemented

    Parameters
    ----------
    weightVal : float, optional
        Set weights to this values (0=flagged), by default 1.
    soltabImport : str, optional
        Name of a soltab. Copy weights from this soltab, by default do not copy.
    flagBad : bool, optional
        Re-apply flags to bad values, by default False.
    """

    import numpy as np

    logging.info("Reweighting soltab: "+soltab.name)

    if mergeSoltab != '':
        solset = soltab.getSolset()
        soltabI = solset.getSoltab(soltabImport)
        soltabI.selection = soltab.selection

        weights, axes = soltab.getValues(weight = True)
        weightsI, axesI = soltabI.getValues(weight = True)
        if axes.keys() != axesI.keys() or weights.shape != weightsI.shape:
            logging.error('Impossible merge two tables with different axes values')
            return 1
        weights[ np.where(weightsI == 0) ] = 0.
        soltab.addHistory('WEIGHT imported from '+soltabI+'.')
    else:
        weights = weightVal
        sw.addHistory('REWEIGHTED to '+str(weightVal)+'.')

    soltab.setValues(weights, weight=True)

    if flagBad:
        weights = soltab.getValues(weight = True, retAxesVals = False)
        vals = soltab.getValues(retAxesVals = False)
        if soltab.getType() == 'amplitude': weights[np.where(vals == 1)] = 0
        else: weights[np.where(vals == 0)] = 0
        soltab.setValues(weights, weight=True)

    return 0
