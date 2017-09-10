#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Residual operation for LoSoTo

# This operation subtract a clock and/or tec from a phase
# Operation is flag-only capable

import logging
from losoto.operations_lib import *

logging.debug('Loading RESIDUALS module.')

def run_parser(soltab, parser, step):
    soltabsToSub = parser.getarraystr( step, 'soltabsToSub' ) # no default
    ratio = parser.getbool( step, 'ratio', False )
    return run(soltab, soltabsToSub, ratio)


def run( soltab, soltabsToSub, ratio=False ):
    """
    Subtract/divide two tables or a clock/tec/tec3rd/rm from a phase.

    Parameters
    ----------
    soltabsToSub : list of str
        List of soltabs to subtract

    ratio : bool, optional
        Return the ratio instead of subtracting, by default False.
    """
    import numpy as np

    logging.info("Subtract soltab: "+soltab.name)

    allsoltabsub = [] # sub tables
    solset = soltab.getSolset()
    for soltabToSub in soltabsToSub:
        soltabsub = solset.getSoltab(soltabToSub)

        # selection
        soltabsub.selection = soltab.selection

        if soltab.getType() != 'phase' and (soltabsub.getType() == 'tec' or soltabsub.getType() == 'clock' or soltabsub.getType() == 'rotationmeasure' or soltabsub.getType() == 'tec3rd'):
            logging.warning(soltabToSub+' is of type clock/tec/rm and should be subtracted from a phase. Skipping it.')
            return 1
        allsoltabsub.append( soltabsub )
        logging.info('Subtracting table: '+soltabToSub)

        # a major speed up if tables are assumed with same axes, check that (should be the case in almost any case)
        for axisName in soltabsub.getAxesNames():
            assert all(soltabsub.getAxisValues(axisName) == soltab.getAxisValues(axisName))
    
    if soltabsub.getType() == 'clock' or soltabsub.getType() == 'tec' or soltabsub.getType() == 'tec3rd' or soltabsub.getType() == 'rotationmeasure':
        # the only return axes is freq, slower but better code
        for vals, weights, coord, selection in soltab.getValuesIter(returnAxes='freq', weight = True):

            for soltabsub in allsoltabsub:

                # restrict to have the same coordinates of phases
                for i, axisName in enumerate(soltabsub.getAxesNames()):
                    soltabsub.selection[i] = selection[soltab.getAxesNames().index(axisName)]

                valsSub = np.squeeze(soltabsub.getValues(retAxesVals=False, weight=False))
                weightsSub = np.squeeze(soltabsub.getValues(retAxesVals=False, weight=True))

                if soltabsub.getType() == 'clock':
                    vals -= 2. * np.pi * valsSub * coord['freq']

                elif soltabsub.getType() == 'tec':
                    vals -= -8.44797245e9 * valsSub / coord['freq']

                elif soltabsub.getType() == 'tec3rd':
                    vals -= - 1.e21 * valsSub / np.power(coord['freq'],3)

                elif soltabsub.getType() == 'rotationmeasure':
                    wav = 2.99792458e8/coord['freq']
                    ph = wav * wav * valsSub
                    if coord['pol'] == 'XX' or coord['pol'] == 'RR':
                        vals -= ph
                    elif coord['pol'] == 'YY' or coord['pol'] == 'LL':
                        vals += ph

                # flag data that are contaminated by flagged clock/tec data
                if weightsSub == 0: weights[:] = 0

            soltab.setValues(vals, selection)
            soltab.setValues(weights, selection, weight = True)
    else:
        for soltabsub in allsoltabsub:
            if ratio: soltab.setValues((soltab.getValues(retAxesVals=False)-soltabsub.getValues(retAxesVals=False))/soltabsub.getValues(retAxesVals=False))
            else: soltab.setValues(soltab.getValues(retAxesVals=False)-soltabsub.getValues(retAxesVals=False))
            weight = soltab.getValues(retAxesVals=False, weight=True)
            weight[soltabsub.getValues(retAxesVals=False, weight=True) == 0] = 0
            soltab.setValues(weight, weight = True)
        
    soltab.addHistory('RESIDUALS by subtracting tables '+' '.join(soltabsToSub))
    soltab.flush()
        
    return 0

