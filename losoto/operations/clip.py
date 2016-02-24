#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Clip solutions around the median by a factor specified by the user.
# WEIGHT: flag compliant, putting weights into median is tricky

# Implemented by Martin Hardcastle based on original flag.py code

import logging
from losoto.operations_lib import *

logging.debug('Loading CLIP module.')

def run( step, parset, H ):

    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    axesToClip = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    clipLevel = parset.getFloat('.'.join(["LoSoTo.Steps", step, "ClipLevel"]), 0. )
    log = parset.getBool('.'.join(["LoSoTo.Steps", step, "Log"]), True )
    
    if len(axesToClip) < 1:
        logging.error("Please specify axes to clip.")
        return 1
    if clipLevel == 0.:
        logging.error("Please specify factor above/below median at which to clip.")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Clipping soltab: "+soltab._v_name)

        sf = solFetcher(soltab)

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        # some checks
        for i, axis in enumerate(axesToClip[:]):
            if axis not in sf.getAxesNames():
                del axesToClip[i]
                logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

        if sf.getType() != 'amplitude':
            logging.error('CLIP is for "amplitude" tables, not %s.' % sf.getType())
            continue

        sw = solWriter(soltab, useCache=True) # remember to flush()

        before_count=0
        after_count=0
        total=0
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=axesToClip, weight = True):

            total+=len(vals)
            before_count+=(len(weights)-np.count_nonzero(weights))

            # first find the median and standard deviation
            if (weights == 0).all():
                valmedian = 0
            else:
                if log:
                    valmedian = np.median(np.log10(vals[(weights != 0)]))
                    rms = np.std(np.log10(vals[(weights != 0)]))
                    np.putmask(weights, np.abs(np.log10(vals)-valmedian) > rms * clipLevel, 0)
                else:
                    valmedian = np.median(vals[(weights != 0)])
                    rms = np.std(vals[(weights != 0)])
                    np.putmask(weights, np.abs(vals-valmedian) > rms * clipLevel, 0)
        
            after_count+=(len(weights)-np.count_nonzero(weights))

            # writing back the solutions
            sw.selection = selection
            sw.setValues(weights, weight=True)

        sw.addHistory('CLIP (over %s with %s sigma cut)' % (axesToClip, clipLevel))
        logging.info('Clip, flagged data: %f %% -> %f %%' \
                % (100.*before_count/total, 100.*after_count/total))

        sw.flush()
        
    return 0


