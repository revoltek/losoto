#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading REFERENCE module.')

def _run_parser(soltab, parser, step):
    refAnt = parser.getstr( step, 'refAnt', '' )
    refDir = parser.getstr( step, 'refDir', '' )

    parser.checkSpelling( step, soltab, ['refAnt','refDir'])
    return run(soltab, refAnt, refDir)


def run( soltab, refAnt='', refDir=''):
    """
    Reference to an antenna

    Parameters
    ----------
    refAnt : str, optional
        Reference antenna for phases or "best" to use the least flagged antenna. Empty string does not change phases. Default: ''.
    refDir : str, optional
        Reference direction for phases. Empty string does not change phases. Default: ''.
    """

    if not soltab.getType() in ['phase', 'scalarphase', 'rotation', 'tec', 'clock', 'tec3rd', 'rotationmeasure']:
        logging.error('Reference possible only for phase, scalarphase, clock, tec, tec3rd, rotation and rotationmeasure solution tables. Ignore referencing.')
        return 1

    if refAnt == 'best': 
        weights = soltab.getValues(retAxesVals=False, weight=True)
        weights = np.sum(weights, axis=tuple([i for i, axis_name in enumerate(soltab.getAxesNames()) if axis_name != 'ant']), dtype=np.float)
        refAnt = soltab.getAxisValues('ant')[np.where(weights == np.max(weights))[0][0]]
        logging.info('Using %s for reference antenna.' % refAnt)

    elif not refAnt in soltab.getAxisValues('ant', ignoreSelection = True) and refAnt != '':
        logging.error('Reference antenna '+refAnt+' not found.')
        return 1

    if not refDir in soltab.getAxisValues('dir', ignoreSelection = True) and refDir != '':
        logging.error('Reference direction '+refDir+' not found.')
        return 1

    # get reference direction if needed
    if refDir != '' and refAnt != '':
        logging.info('Referencing on both dir and ant.')
        soltab.setSelection(dir=refDir, ant=refAnt)
        valsRef = soltab.getValues(retAxesVals=False)
        soltab.clearSelection()
        vals = soltab.getValues(retAxesVals=False)
        dirAxis = soltab.getAxesNames().index('dir')
        antAxis = soltab.getAxesNames().index('ant')
        print ('len vals ref', valsRef.shape)
        valsRef = np.repeat(valsRef, axis=dirAxis, repeats=len(soltab.getAxisValues('dir')))
        valsRef = np.repeat(valsRef, axis=antAxis, repeats=len(soltab.getAxisValues('ant')))
        vals = vals - valsRef

    elif refDir != '':
        soltab.setSelection(dir=refDir)
        valsRefDir = soltab.getValues(retAxesVals=False)
        soltab.clearSelection()
        vals = soltab.getValues(retAxesVals=False)
        dirAxis = soltab.getAxesNames().index('dir')
        vals = vals - np.repeat(valsRefDir, axis=dirAxis, repeats=len(soltab.getAxisValues('dir')))

    # use automatic antenna referencing
    elif refAnt != '':
        vals = soltab.getValues(retAxesVals=False, refAnt=refAnt)

    soltab.setValues(vals)

    soltab.addHistory('REFERENCED (to antenna: %s)' % (refAnt))
    return 0
