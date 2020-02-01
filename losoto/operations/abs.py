#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading ABS module.')

def _run_parser(soltab, parser, step):
    parser.checkSpelling( step, soltab )
    return run(soltab)

def run( soltab ):
    """
    Take absolute value. Needed before smooth if amplitudes are negative!
    WEIGHT: no need to be weight compliant

    Parameters
    ----------
    soltab : soltab obj
        Solution table.
    """

    import numpy as np

    logging.info("Taking ABSolute value of soltab: "+soltab.name)

    vals = soltab.getValues(retAxesVals = False)
    count = np.count_nonzero(vals<0)

    logging.info('Abs: %i points initially negative (%f %%)' % (count,100*float(count)/np.count_nonzero(vals)))

    # writing back the solutions
    soltab.setValues(np.abs(vals))

    soltab.addHistory('ABSolute value taken')
        
    return 0


