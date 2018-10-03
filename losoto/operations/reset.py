#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
import logging

logging.debug('Loading RESET module.')

def _run_parser(soltab, parser, step):
    dataVal = parser.getfloat( step, 'dataVal' ) # no default

    return run(soltab, dataVal=dataVal)

def run( soltab, dataVal=None ):
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

    if dataVal is None:
        if solType == 'amplitude':
            dataVal = 1.
        else:
            dataVal = 0.
    
    soltab.setValues(dataVal)

    soltab.addHistory('RESET')
    return 0
