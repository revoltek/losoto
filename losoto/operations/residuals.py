#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Residual operation for LoSoTo

# This operation subtract a clock and/or tec from a phase
# Operation is flag-only capable

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading RESIDUALS module.')

def _run_parser(soltab, parser, step):
    soltabsToSub = parser.getarraystr( step, 'soltabsToSub' ) # no default
    ratio = parser.getbool( step, 'ratio', False )

    parser.checkSpelling( step, soltab, ['soltabsToSub','ratio'])
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

    solset = soltab.getSolset()
    for soltabToSub in soltabsToSub:
        soltabsub = solset.getSoltab(soltabToSub)

        if soltab.getType() != 'phase' and (soltabsub.getType() == 'tec' or soltabsub.getType() == 'clock' or soltabsub.getType() == 'rotationmeasure' or soltabsub.getType() == 'tec3rd'):
            logging.warning(soltabToSub+' is of type clock/tec/rm and should be subtracted from a phase. Skipping it.')
            return 1
        logging.info('Subtracting table: '+soltabToSub)

        # a major speed up if tables are assumed with same axes, check that (should be the case in almost any case)
        for i, axisName in enumerate(soltabsub.getAxesNames()):
            # also armonise selection by copying only the axes present in the outtable and in the right order
            soltabsub.selection[i] = soltab.selection[soltab.getAxesNames().index(axisName)]
            assert (soltabsub.getAxisValues(axisName) == soltab.getAxisValues(axisName)).all() # table not conform

        if soltab.getValues(retAxesVals=False, weight=False).shape != soltabsub.getValues(retAxesVals=False, weight=False).shape:
            hasMissingAxes = True
        else:
            hasMissingAxes = False

        if soltabsub.getType() == 'clock' or soltabsub.getType() == 'tec' or soltabsub.getType() == 'tec3rd' or soltabsub.getType() == 'rotationmeasure' or hasMissingAxes:

            freq = soltab.getAxisValues('freq')
            vals = soltab.getValues(retAxesVals=False, weight=False)
            weights = soltab.getValues(retAxesVals=False, weight=True)
            #print 'vals', vals.shape

            # valsSub doesn't have freq
            valsSub = soltabsub.getValues(retAxesVals=False, weight=False)
            weightsSub = soltabsub.getValues(retAxesVals=False, weight=True)
            #print 'valsSub', valsSub.shape

            # add missing axes and move it to the last position
            expand = [soltab.getAxisLen(ax) for ax in soltab.getAxesNames() if ax not in soltabsub.getAxesNames()]
            #print "expand:", expand
            valsSub = np.resize( valsSub, expand+list(valsSub.shape) )
            weightsSub = np.resize( weightsSub, expand+list(weightsSub.shape) )
            #print 'valsSub missing axes', valsSub.shape

            # reorder axes to match soltab
            names = [ax for ax in soltab.getAxesNames() if ax not in soltabsub.getAxesNames()] + soltabsub.getAxesNames()
            #print names, soltab.getAxesNames()
            valsSub = reorderAxes( valsSub, names, soltab.getAxesNames() )
            weightsSub = reorderAxes( weightsSub, names, soltab.getAxesNames() )
            weights[ weightsSub == 0 ] = 0 # propagate flags
            #print 'valsSub reorder', valsSub.shape

            # put freq axis at the end
            idxFreq = soltab.getAxesNames().index('freq')
            vals = np.swapaxes(vals, idxFreq, len(vals.shape)-1)
            valsSub = np.swapaxes(valsSub, idxFreq, len(valsSub.shape)-1)
            #print 'vals reshaped', valsSub.shape

            # a multiplication will go along the last axis of the array
            if soltabsub.getType() == 'clock':
                vals -= 2. * np.pi * valsSub * freq

            elif soltabsub.getType() == 'tec':
                vals -= -8.44797245e9 * valsSub / freq

            elif soltabsub.getType() == 'tec3rd':
                vals -= - 1.e21 * valsSub / np.power(freq,3)

            elif soltabsub.getType() == 'rotationmeasure':
                # put pol axis at the beginning
                idxPol = soltab.getAxesNames().index('pol')
                if idxPol == len(vals.shape)-1: idxPol = idxFreq # maybe freq swapped with pol
                vals = np.swapaxes(vals, idxPol, 0)
                valsSub = np.swapaxes(valsSub, idxPol, 0)
                #print 'vals reshaped 2', valsSub.shape

                wav = 2.99792458e8/freq
                ph = wav * wav * valsSub
                #if coord['pol'] == 'XX' or coord['pol'] == 'RR':
                pols = soltab.getAxisValues('pol')
                assert len(pols) == 2 # full jons not supported
                if (pols[0] == 'XX' and pols[1] == 'YY') or \
                   (pols[0] == 'RR' and pols[1] == 'LL'):
                        vals[0] -= ph[0]
                        vals[1] += ph[1]
                else:
                        vals[0] += ph[0]
                        vals[1] -= ph[1]

                vals = np.swapaxes(vals, 0, idxPol)

            else:
                if ratio:
                    vals = (vals - valsSub) / valsSub
                else:
                    vals -= valsSub

            # move freq axis back
            vals = np.swapaxes(vals, len(vals.shape)-1, idxFreq)

            soltab.setValues(vals)
            soltab.setValues(weights, weight = True)
        else:
            if ratio: soltab.setValues((soltab.getValues(retAxesVals=False)-soltabsub.getValues(retAxesVals=False))/soltabsub.getValues(retAxesVals=False))
            else: soltab.setValues(soltab.getValues(retAxesVals=False)-soltabsub.getValues(retAxesVals=False))
            weight = soltab.getValues(retAxesVals=False, weight=True)
            weight[soltabsub.getValues(retAxesVals=False, weight=True) == 0] = 0
            soltab.setValues(weight, weight = True)

    soltab.addHistory('RESIDUALS by subtracting tables '+' '.join(soltabsToSub))

    return 0

