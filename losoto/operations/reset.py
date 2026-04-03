#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading RESET module.')

def _run_parser(soltab, parser, step):
    dataVal = parser.getfloat( step, 'dataVal', -999. )

    parser.checkSpelling( step, soltab, ['dataVal'])
    return run(soltab, dataVal)

def run( soltab, dataVal=-999. ):
    """
    This operation reset all the selected solution values.
    WEIGHT: flag compliant, no need for weight

    Parameters
    ----------
    dataVal : float, optional
        If given set values to this number, otherwise uses 1 for amplitude and 0 for all other soltab types.

    Example : to reset a fulljones amplitude
              [resetA]
              operation=RESET
              soltab=sol000/amplitude000
              pol=XX,YY
              dataVal=1.
              [resetB]
              operation=RESET
              soltab=sol000/amplitude000
              pol=XY,YX
              dataVal=0.
    """

    logging.info("Resetting soltab: "+soltab.name)

    solType = soltab.getType()

    if dataVal == -999.:
        if solType == 'amplitude':
            axesNames = soltab.getAxesNames()
            # special case full-Jones matrix
            if 'pol' in axesNames:
                polAxis = axesNames.index('pol')
                if soltab.val.shape[polAxis] == 4: # assume four dimensions == full Jones
                    logging.info("Detected full-Jones amplitudes")
                    val, axesval = soltab.getValues()
                    pols = axesval['pol']
                    # iterate polarizations and set to one/zero depending on diagonal/off-diagonal
                    for i, pol in enumerate(pols):
                        # Build a slicing tuple
                        slicer = [slice(None)] * val.ndim
                        slicer[polAxis] = i # in the polarization dimension, pick the i-th entry
                        resetval = 1. if pol in ['XX', 'YY', 'RR', 'LL'] else 0.
                        # logging.info(f"Reset pol {pol} to {resetval}")
                        val[tuple(slicer)] = resetval
                    print(val.shape)
                    soltab.setValues(val)
                    soltab.addHistory('RESET')
                    return 0

            dataVal = 1.
        else:
            dataVal = 0.
    
    soltab.setValues(dataVal)

    soltab.addHistory('RESET')
    return 0
