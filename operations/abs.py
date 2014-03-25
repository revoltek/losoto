#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Take absolute value. Needed before smooth if amplitudes are negative!

# Implemented by Martin Hardcastle based on clip/flag code

import logging
from operations_lib import *

logging.debug('Loading ABS module.')

def run( step, parset, H ):

    import numpy as np
    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    # No need to specify an axis, just use time
    axesToAbs = ['time']

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Taking ABSolute value of soltab: "+soltab._v_name)

        sf.setSelection(ant=ants, pol=pols, dir=dirs)

        total=0
        count=0
        for vals, coord in sf.getValuesIter(returnAxes=axesToAbs):

            total+=len(vals)
            count+=np.count_nonzero(vals<0)
            valsnew=np.abs(vals)

            # writing back the solutions
            coord = removeKeys(coord, axesToAbs)
            sw.setSelection(**coord)
            sw.setValues(valsnew)

        logging.info('Abs: %i points initially negative (%f %%)' % (count,float(count)/total))

        sw.addHistory('ABSolute value taken')
        
    return 0


