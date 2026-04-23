#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading REFERENCE module.')

def _run_parser(soltab, parser, step):
    refAnt = parser.getstr( step, 'refAnt', '' )
    refDir = parser.getstr( step, 'refDir', '' )
    preservePol = parser.getbool( step, 'preservePol', False)

    parser.checkSpelling( step, soltab, ['refAnt','refDir','preservePol'])
    return run(soltab, refAnt, refDir, preservePol)


def run( soltab, refAnt='', refDir='', preservePol=False):
    """
    Reference to an antenna

    Parameters
    ----------
    refAnt : str, optional
        Reference antenna for phases. Empty string does not change phases. Default: ''.
    refDir : str, optional
        Reference direction for phases. Empty string does not change phases. Default: ''.
    preservePol : bool, optional
        For fullJones, or if preservePol is given, reference only to X/R polarization to preserve global X-Y delay
    """

    if not soltab.getType() in ['phase', 'scalarphase', 'rotation', 'tec', 'clock', 'tec3rd', 'rotationmeasure']:
        logging.error('Reference possible only for phase, scalarphase, clock, tec, tec3rd, rotation and rotationmeasure solution tables. Ignore referencing.')
        return 1

    if refAnt != 'closest' and refAnt != 'auto' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True) and refAnt != '':
        logging.error('Reference antenna '+refAnt+' not found.')
        return 1

    if refDir != '':
        if not refDir in soltab.getAxisValues('dir', ignoreSelection = True):
            logging.error('Reference direction '+refDir+' not found.')
            return 1

    if refDir != '' and refAnt != '':
        vals = soltab.getValues(retAxesVals=False, refAnt=refAnt, refDir=refDir, preservePol=preservePol)

    elif refDir != '':
        vals = soltab.getValues(retAxesVals=False, refDir=refDir, preservePol=preservePol)

    elif refAnt != '':
        vals = soltab.getValues(retAxesVals=False, refAnt=refAnt, preservePol=preservePol)

    soltab.setValues(vals)

    if refAnt != '': soltab.addHistory('REFERENCED (to antenna: %s, preservePol: %s)' % (refAnt, preservePol))
    if refDir != '': soltab.addHistory('REFERENCED (to direction: %s)' % (refDir))
    return 0
