#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a smoothing function

import logging
from operations_lib import *

logging.debug('Loading SMOOTH module.')

def run( step, parset, H ):

    import scipy.ndimage.filters
    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    axesToSmooth = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    FWHM = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "FWHM"]), [] )

    if len(axesToSmooth) != len(FWHM):
        logging.error("Axes and FWHM lenghts must be equal.")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Smoothing soltab: "+soltab._v_name)

        sf.setSelection(ant=ants, pol=pols, dir=dirs)

        for i, axis in enumerate(axesToSmooth[:]):
            if axis not in sf.getAxesNames():
                del axesToSmooth[i]
                del FWHM[i]
                logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

        for vals, coord in sf.getValuesIter(returnAxes=axesToSmooth):

            valsnew = scipy.ndimage.filters.median_filter(vals, FWHM)
        
            # TODO: implement flag control
            # find the local mean (without any nans)
#            valmean = np.ma.mean(np.ma.masked_array(vals, np.isnan(vals)), axis=0)
            # print 'mean value: ' + str(valmean)
            # replace any nans with median
#            valbool = np.isnan(vals)
#            nanindex = np.where(valbool)
#            if valbool.any():
#                logging.warning('Found NaNs in solutions. Replacing with local mean.')
                # check if valmean is iterable
#                if valmean.shape is tuple():
#                    np.putmask(vals, valbool, valmean)
#                else:
#                    for x,mval in enumerate(valmean):
#                        np.putmask(vals[:,x], valbool[:,x], mval)

            # writing back the solutions (selection on all the coord axis)
            # this can be properly broacasted
            coord = removeKeys(coord, axesToSmooth)
            sw.setSelection(**coord)
            sw.setValues(valsnew)

        sw.addHistory('SMOOTH (over %s with box size = %s)' % (axesToSmooth, FWHM))
    return 0


