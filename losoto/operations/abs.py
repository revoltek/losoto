#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Take absolute value. Needed before smooth if amplitudes are negative!
# WEIGHT: no need to be weight compliant

# Implemented by Martin Hardcastle based on clip/flag code

import logging
from losoto.operations_lib import *

logging.debug('Loading ABS module.')

def run( step, parset, H ):

    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Taking ABSolute value of soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        vals = sf.getValues(retAxesVals = False)
        count = np.count_nonzero(vals<0)

        logging.info('Abs: %i points initially negative (%f %%)' % (count,100*float(count)/np.count_nonzero(vals)))

        # writing back the solutions
        sw.setValues(np.abs(vals))

        sw.addHistory('ABSolute value taken')
        
    return 0


