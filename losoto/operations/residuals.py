#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Residual operation for LoSoTo

# This operation subtract a clock and/or tec from a phase
# Operation is flag-only capable

import logging
from losoto.operations_lib import *

logging.debug('Loading RESIDUALS module.')

def run( step, parset, H ):
    """
    subtract two tables or a clock/tec/tec3rd/rm from a phase.
    """
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    soltabsToSub = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Sub"]), [] )
    ratio = parset.getBool('.'.join(["LoSoTo.Steps", step, "Ratio"]), False )

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab, useCache = True)

        # this will make a selection for the getValues() and getValuesIter()
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)
        sw.setSelection(**userSel)

        sfss = [] # sol fetcher to sub tables
        for soltabToSub in soltabsToSub:
            ss, st = soltabToSub.split('/')
            sfs = solFetcher(H.getSoltab(ss, st))

            # selection
            userSel = {}
            for axis in sfs.getAxesNames():
                userSel[axis] = getParAxis( step, parset, H, axis )
            sfs.setSelection(**userSel)

            if sf.getType() != 'phase' and (sfs.getType() == 'tec' or sfs.getType() == 'clock' or sfs.getType() == 'rotationmeasure' or sfs.getType() == 'tec3rd'):
                logging.warning(soltabToSub+' is of type clock/tec/rm and should be subtracted from a phase. Skipping it.')
                continue
            sfss.append( sfs )
            logging.info('Subtracting table: '+soltabToSub)

            # a major speed up if tables are assumed with same axes, check that (should be the case in almost any case)
            for axisName in sfs.getAxesNames():
                assert all(sfs.getAxisValues(axisName) == sf.getAxisValues(axisName))
        
        if sf.getType() == 'phase' and (sfs.getType() == 'tec' or sfs.getType() == 'clock' or sfs.getType() == 'rotationmeasure' or sfs.getType() == 'tec3rd' ):
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

                    elif sfs.getType() == 'tec3rd':
                        vals -= - 1.e21 * valsSub / np.power(coord['freq'],3)

                    elif sfs.getType() == 'rotationmeasure':
                        wav = 2.99792458e8/coord['freq']
                        ph = wav * wav * valsSub
                        if coord['pol'] == 'XX' or coord['pol'] == 'RR':
                            vals -= ph
                        elif coord['pol'] == 'YY' or coord['pol'] == 'LL':
                            vals += ph
                    else:
                        vals -= valsSub

                    # flag data that are contaminated by flagged clock/tec data
                    if weightsSub == 0: weights[:] = 0

                sw.selection = selection
                sw.setValues(vals)
                sw.setValues(weights, weight = True)
        else:
                if ratio: sw.setValues((sf.getValues(retAxesVals=False)-sfs.getValues(retAxesVals=False))/sfs.getValues(retAxesVals=False))
                else: sw.setValues(sf.getValues(retAxesVals=False)-sfs.getValues(retAxesVals=False))
                weight = sf.getValues(retAxesVals=False, weight=True)
                weight[sfs.getValues(retAxesVals=False, weight=True) == 0] = 0
                sw.setValues(weight, weight = True)
            
        sw.addHistory('RESIDUALS by subtracting tables '+' '.join(soltabsToSub))
        sw.flush()
        del sf
        del sw
        
    return 0

