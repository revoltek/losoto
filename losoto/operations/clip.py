#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

logging.debug('Loading CLIP module.')

def _run_parser(soltab, parser, step):
    axesToClip = parser.getarraystr( step, 'axesToClip' ) # no default
    clipLevel = parser.getfloat( step, 'clipLevel', 5. )
    log = parser.getbool( step, 'log', True )
    return run(soltab, axesToClip, clipLevel, log)

def run( soltab, axesToClip, clipLevel=5., log=True ):
    """
    Clip solutions around the median by a factor specified by the user.
    WEIGHT: flag compliant, putting weights into median is tricky

    Parameters
    ----------
    axesToClip : list of str
        axes along which to calculate the median

    clipLevel : float, optional
        factor above/below median at which to clip, by default 5

    log : bool, optional
        clip is done in log10 space, by default True
    """

    import numpy as np

    logging.info("Clipping soltab: "+soltab.name)

    # input check
    if len(axesToClip) < 1:
        logging.error("Please specify axes to clip.")
        return 1

    if clipLevel <= 0.:
        logging.error("Please specify a positive factor above/below median at which to clip.")
        return 1

    # some checks
    for i, axis in enumerate(axesToClip[:]):
        if axis not in soltab.getAxesNames():
            del axesToClip[i]
            logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

    if soltab.getType() != 'amplitude':
        logging.error('CLIP is for "amplitude" tables, not %s.' % soltab.getType())
        return 1

    before_count=0
    after_count=0
    total=0
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToClip, weight = True):

        total += len(vals)
        before_count += (len(weights)-np.count_nonzero(weights))

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
    
        after_count += (len(weights)-np.count_nonzero(weights))

        # writing back the solutions
        soltab.setValues(weights, selection, weight=True)

    soltab.addHistory('CLIP (over %s with %s sigma cut)' % (axesToClip, clipLevel))
    logging.info('Clip, flagged data: %f %% -> %f %%' \
            % (100.*before_count/total, 100.*after_count/total))

    soltab.flush()
        
    return 0


