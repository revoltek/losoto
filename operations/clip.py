#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Clip solutions around the median by a factor specified by the user.
# Flagging is implemented by setting values to NaN.

# Implemented by Martin Hardcastle based on original flag.py code

import logging
from operations_lib import *

logging.debug('Loading CLIP module.')

def run( step, parset, H ):

    import numpy as np
    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    axesToClip = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    clipLevel = parset.getFloat('.'.join(["LoSoTo.Steps", step, "ClipLevel"]), 0. )
    
    if len(axesToClip) < 1:
        logging.error("Please specify axes to clip.")
        return 1
    if clipLevel == 0.:
        logging.error("Please specify factor above/below median at which to clip.")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Clipping soltab: "+soltab._v_name)

        sf.setSelection(ant=ants, pol=pols, dir=dirs)

        # some checks

        for i, axis in enumerate(axesToClip[:]):
            if axis not in sf.getAxesNames():
                del axesToClip[i]
                logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

        before_count=0
        after_count=0
        total=0
        for vals, coord in sf.getValuesIter(returnAxes=axesToClip):

            total+=len(vals)
            before_count+=np.count_nonzero(np.isnan(vals))

            # clipping
            # first find the median and standard deviation
            valmedian = np.median(vals)
            clipvalue = valmedian * clipLevel
            np.putmask(vals, vals > clipvalue, np.nan)
            clipvalue = valmedian / clipLevel
            np.putmask(vals, vals < clipvalue, np.nan)
        
            after_count+=np.count_nonzero(np.isnan(vals))

            # writing back the solutions
            coord = removeKeys(coord, axesToClip)
            sw.setSelection(**coord)
            sw.setValues(vals)

        sw.addHistory('CLIP (over %s with %s sigma cut)' % (axesToClip, clipLevel))
        logging.info('Clip: %i points initially bad, %i after clipping (%f %%)' % (before_count,after_count,float(after_count)/total))
        
    return 0


