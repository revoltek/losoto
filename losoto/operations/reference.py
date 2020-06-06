#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading REFERENCE module.')

def _run_parser(soltab, parser, step):
    refAnt = parser.getstr( step, 'refAnt', '' )

    parser.checkSpelling( step, soltab, ['refAnt'])
    return run(soltab, refAnt)


def run( soltab, refAnt=''):
    """
    Reference to an antenna

    Parameters
    ----------
    refAnt : str, optional
        Reference antenna for phases or "best" to use the least flagged antenna. Default: best.
    """

    if not soltab.getType() in ['phase', 'scalarphase', 'rotation', 'tec', 'clock', 'tec3rd', 'rotationmeasure']:
        logging.error('Reference possible only for phase, scalarphase, clock, tec, tec3rd, rotation and rotationmeasure solution tables. Ignore referencing.')
        return 1

    if refAnt == '' or refAnt == 'best': 
        weights = soltab.getValues(retAxesVals=False, weight=True)
        weights = np.sum(weights, axis=tuple([i for i, axis_name in enumerate(soltab.getAxesNames()) if axis_name != 'ant']), dtype=np.float)
        refAnt = soltab.getAxisValues('ant')[np.where(weights == np.max(weights))[0][0]]
        logging.info('Using %s for reference antenna.' % refAnt)

    elif not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.error('Reference antenna '+refAnt+' not found.')
        return 1

    vals = soltab.getValues(retAxesVals=False, reference=refAnt)
    soltab.setValues(vals)

    soltab.addHistory('REFERENCED (to antenna: %s)' % (refAnt))
    return 0
