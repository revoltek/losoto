#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
import logging

logging.debug('Loading RESET module.')

def _run_parser(soltab, parser, step):
    return run(soltab)

def run( soltab ):
    """
    This operation reset all the selected amplitudes to 1
    and all other selected solution types to 0
    WEIGHT: flag compliant, no need for weight
    """

    logging.info("Resetting soltab: "+soltab.name)

    solType = soltab.getType()

    if solType == 'amplitude':
        soltab.setValues(1.)
    else:
        soltab.setValues(0.)

    soltab.addHistory('RESET')
    return 0
