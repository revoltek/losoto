#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Residual operation for LoSoTo

# This operation subtract a clock and/or tec from a phase
# Operation is flag-only capable

import logging
from losoto.operations_lib import *

logging.debug('Loading RESIDUAL module.')

def run( step, parset, H ):
    """
    subtract a clock and/or tec from a phase.
    """
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    soltabsToSub = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Sub"]), [] )

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab, useCache = True)

        if sf.getType() != 'phase':
            logging.warning(soltab._v_name+' is not of type phase, ignore.')
            continue

        sfss = [] # sol fetcher to sub tables
        for soltabToSub in soltabsToSub:
            ss, st = soltabToSub.split('/')
            sfs = solFetcher(H.getSoltab(ss, st))
            if sfs.getType() != 'tec' and sfs.getType() != 'clock' and sfs.getType() != 'rotationmeasure':
                logging.warning(soltabToSub+' is not of type clock/tec/rm and cannot be subtracted, ignore.')
                continue
            sfss.append( sfs )
            logging.info('Subtracting table: '+soltabToSub)

            # a major speed up if tables are assumed with same axes, check that (should be the case in almost any case)
            for axisName in sfs.getAxesNames():
                assert all(sfs.getAxisValues(axisName) == sf.getAxisValues(axisName))
        
        # the only return axes is freq, slower but better code
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes='freq', weight = True):

            for sfs in sfss:

                # restrict to have the same coordinates of phases
                for i, axisName in enumerate(sfs.getAxesNames()):
                    sfs.selection[i] = selection[sf.getAxesNames().index(axisName)]

                valsSub = np.squeeze(sfs.getValues(retAxesVals=False, weight=False))
                weightsSub = np.squeeze(sfs.getValues(retAxesVals=False, weight=True))

                if sfs.getType() == 'clock':
                    vals -= 2. * np.pi * valsSub * coord['freq']

                elif sfs.getType() == 'tec':
                    vals -= -8.44797245e9 * valsSub / coord['freq']

                elif sfs.getType() == 'rotationmeasure':
                    wav = 2.99792458e8/coord['freq']
                    ph = wav * wav * valsSub
                    if coord['pol'] == 'XX':
                        vals -= ph
                    elif coord['pol'] == 'YY':
                        vals += ph

                # flag data that are contaminated by flagged clock/tec data
                weights[np.where(weightsSub == 0)] == 0

            sw.selection = selection
            sw.setValues(vals)
            sw.setValues(weights, weight = True)
            
        sw.addHistory('RESIDUALS by subtracting tables '+' '.join(soltabsToSub))
        sw.flush()
        del sf
        del sw
        
    return 0

