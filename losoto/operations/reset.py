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
    """

    logging.info("Resetting soltab: "+soltab.name)

    solType = soltab.getType()

    if dataVal == -999.:
        if solType == 'amplitude':
            dataVal = 1.
        else:
            dataVal = 0.
    
    soltab.setValues(dataVal)

    soltab.addHistory('RESET')
    return 0
