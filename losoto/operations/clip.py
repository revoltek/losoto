#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading CLIP module.')

def _run_parser(soltab, parser, step):
    axesToClip = parser.getarraystr( step, 'axesToClip', [] )
    clipLevel = parser.getfloat( step, 'clipLevel', 5. )
    log = parser.getbool( step, 'log', False )
    mode = parser.getstr( step, 'mode', 'median' )
    
    parser.checkSpelling( step, soltab, ['axesToClip', 'clipLevel', 'log', 'mode'] )
    return run(soltab, axesToClip, clipLevel, log, mode)

def run( soltab, axesToClip=None, clipLevel=5., log=False, mode='median' ):
    """
    Clip solutions around the median by a factor specified by the user.
    WEIGHT: flag compliant, putting weights into median is tricky

    Parameters
    ----------
    axesToClip : list of str
        axes along which to calculate the median (e.g. [time,freq]).

    clipLevel : float, optional
        factor above/below median at which to clip, by default 5.

    log : bool, optional
        clip is done in log10 space, by default False.

    mode : str, optional
        if "median" then flag at rms*clipLevel times away from the median, if it is "above", 
        then flag all values above clipLevel, if it is "below" then flag all values below clipLevel. 
        By default median.
    """

    import numpy as np

    def percentFlagged(w):
        return 100.*(weights.size-np.count_nonzero(weights))/float(weights.size)

    logging.info("Clipping soltab: "+soltab.name)

    # input check
    if len(axesToClip) < 1 and mode == 'median':
        logging.error("Please specify axes to clip.")
        return 1
    else:
        axesToClip = soltab.getAxesNames()

    if mode != 'median' and mode != 'above' and mode != 'below':
        logging.error("Mode can be only: median, above or below.")
        return 1

    if mode == 'median' and soltab.getType() == 'amplitude' and not log:
        logging.warning('Amplitude solution tab detected and log=False. Amplitude solution tables should be treated in log space.')

    # some checks
    for i, axis in enumerate(axesToClip[:]):
        if axis not in soltab.getAxesNames():
            del axesToClip[i]
            logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToClip, weight = True):

        initPercent = percentFlagged(weights)

        # skip all flagged
        if (weights == 0).all():
            continue
        
        if log: vals = np.log10(vals)

        if mode == 'median':
            valmedian = np.nanmedian(vals[(weights != 0)])
            rms = np.nanstd(vals[(weights != 0)])
            np.putmask(weights, np.abs(vals-valmedian) > rms * clipLevel, 0)

        elif mode == 'above':
            np.putmask(weights, vals > clipLevel, 0)

        elif mode == 'below':
            np.putmask(weights, vals < clipLevel, 0)
    
        # writing back the solutions
        if log: vals = 10**vals
        soltab.setValues(weights, selection, weight=True)

        #print('max', np.max(vals[(weights != 0)]))
        #print('median', np.nanmedian(vals[(weights != 0)]))
        logging.debug('Percentage of data flagged (%s): %.3f%% -> %.3f%%' \
            % (removeKeys(coord, axesToClip), initPercent, percentFlagged(weights)))

    soltab.addHistory('CLIP (over %s with %s sigma cut)' % (axesToClip, clipLevel))

    soltab.flush()
        
    return 0


