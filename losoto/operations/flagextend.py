#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing

from losoto._logging import logger as logging
from losoto.lib_operations import nproc, removeKeys

logging.debug('Loading FLAGEXTEND module.')

def _run_parser(soltab, parser, step):
    axesToExt = parser.getarraystr( step, 'axesToExt') # no default
    size = parser.getarrayint( step, 'size' ) # no default
    percent = parser.getfloat( step, 'percent', 50. )
    maxCycles = parser.getint( step, 'maxCycles', 3 )
    ncpu = parser.getint( '_global', 'ncpu', 0 )

    parser.checkSpelling( step, soltab, ['axesToExt', 'size', 'percent', 'maxCycles'])
    return run(soltab, axesToExt, size, percent, maxCycles, ncpu)


def _flag(weights, coord, axesToExt, selection, percent=50, size=[0], maxCycles=3):
        """
        Flag data if surreounded by other flagged data
        weights = the weights to convert into flags
        percent = percent of surrounding flagged point to extend the flag

        return: flags array and final rms
        """
        def extendFlag(flags, percent):
            #flags = flags.astype(int)
            if float(np.sum( flags ))/len(flags) > percent/100.:
                return 1
            else:
                return 0

        def percentFlagged(w):
            return 100.*(weights.size-np.count_nonzero(weights))/float(weights.size)


        import numpy as np
        import scipy.ndimage
        initPercent = percentFlagged(weights)

        # if size=0 then extend to all 2*axis, this otherwise create issues with mirroring
        for i, s in enumerate(size):
            if s == 0: size[i] = 2*weights.shape[i]

        for cycle in range(maxCycles):
            flag = scipy.ndimage.filters.generic_filter((weights==0), extendFlag, size=size, mode='mirror', cval=0.0, origin=0, extra_keywords={'percent':percent})
            weights[ ( flag == 1 ) ] = 0
            # no new flags
            if cycle != 0 and np.count_nonzero(flag) == oldFlagCount: break
            oldFlagCount = np.count_nonzero(flag)

        if percentFlagged(weights) == initPercent:
            logging.debug('Percentage of data flagged (%s): %.3f -> None' \
                    % (removeKeys(coord, axesToExt), initPercent))
        else:
            logging.debug('Percentage of data flagged (%s): %.3f -> %.3f %%' \
                    % (removeKeys(coord, axesToExt), initPercent, percentFlagged(weights)))

        return weights, selection


def run( soltab, axesToExt, size, percent=50., maxCycles=3, ncpu=0 ):
    """
    This operation for LoSoTo implement a extend flag procedure
    It can work in multi dimensional space and for each datum check if the surrounding data are flagged to a certain %, then flag also that datum
    The size of the surrounding footprint can be tuned
    WEIGHT: compliant

    Parameters
    ----------
    axesToExt : list of str
        Axes used to find close flags.

    size : list of int
        Size of the window (diameter), per axis. If 0 is given then the entire length of the axis is assumed.
        Must be a vector of same length of Axes.

    percent : float, optional
        Percent of flagged data around the point to flag it, by default 50.

    maxCycles : int, optional
        Number of independent cycles of flag expansion, by default 3.

    ncpu : int, optional
        Number of CPU used, by default all available.
    """

    import numpy as np

    logging.info("Extending flag on soltab: "+soltab.name)

    # input check
    if axesToExt == []:
        logging.error("Please specify at least one axis to extend flag.")
        return 1

    for axisToExt in axesToExt:
        if axisToExt not in soltab.getAxesNames():
            logging.error('Axis \"'+axisToExt+'\" not found.')
            return 1

    args = [
        (weights, coord, axesToExt, selection, percent, size, maxCycles)
        for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=axesToExt, weight=True)
    ]

    ncpu = ncpu if ncpu > 0 else nproc()  # use all available CPUs if ncpu is not set
    logging.debug("Using %s CPU(s) for operation FLAGEXTEND.", ncpu)
    with multiprocessing.Pool(ncpu) as pool:
        results = pool.starmap(_flag, args)

    logging.info('Writing solutions')
    for w, sel in results:
        soltab.setValues(w, sel, weight=True)

    soltab.addHistory('FLAG EXTENDED (over %s)' % (str(axesToExt)))
    return 0
