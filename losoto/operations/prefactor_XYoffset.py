#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

logging.debug('Loading PREFACTOR_XYOFFSET module.')

def _run_parser(soltab, parser, step):
    return run(soltab)


def normalize(phase):
    """
    Normalize phase to the range [-pi, pi].
    """

    # Convert to range [-2*pi, 2*pi].
    out = np.fmod(phase, 2.0 * np.pi)

    # Convert to range [-pi, pi]
    out[out < -np.pi] += 2.0 * np.pi
    out[out > np.pi] -= 2.0 * np.pi

    return out


def run( soltab):

    import numpy as np
    import scipy.signal as s
    import pylab as pl

    logging.info("Running XYoffset on: "+soltab.name)
    solset = soltab.getSolset()
    
    solType = soltab.getType()
    if solType != 'phase':
       logging.warning("Soltab type of "+soltab.name+" is: "+solType+" should be phase. Ignoring.")
       return 1
    
    phases_tmp = np.copy(soltab.val)
    freqs = np.copy(soltab.freq)
    npol = len(soltab.pol)
    sourceID = 0
    refstationID=2
    
    stationsnames = [ stat for stat in soltab.ant]
    # this gets the subband number to any given frequency in HBA-low
    subbands = np.unique(np.round(freqs/195.3125e3-512.))
    nsubbands = len(subbands)
    nchan = len(freqs)/nsubbands
    if nsubbands*nchan != len(freqs):
        print "find_cal_global_phaseoffset_losoto.py: irregular number of ch/SB detected! Bailing out!"
        print "  nchan %d, nSB: %d, nfreq: %d" % (nchan, nsubbands, len(freqs))
        sys.exit(1)
    tmpfreqs = freqs.reshape([nsubbands,nchan])
    freq_per_sb = np.mean(tmpfreqs,axis=1)
    nstations = len(stationsnames)
    refphases = phases_tmp[:,sourceID,refstationID,:,:]

    for istat in xrange(nstations):
        phases_00 = phases_tmp[0,sourceID,istat,:,:]-refphases[0,:,:]
        phases_11 = phases_tmp[1,sourceID,istat,:,:]-refphases[1,:,:]
        phases_diff = normalize(phases_00-phases_11)
        tmp_phases_diff = np.median(phases_diff,axis=1)
        med_phases_diff = np.median(tmp_phases_diff.reshape([nsubbands,nchan]),axis=1)
        if istat == 0:
            global_stat_offsets = med_phases_diff
        else:
            global_stat_offsets = np.vstack( (global_stat_offsets, med_phases_diff) )
    global_stat_offsets_smoothed = np.zeros([nsubbands,nstations,npol])
    for istat in xrange(nstations):
        global_stat_offsets_smoothed[:,istat,-1] = s.medfilt(global_stat_offsets[istat,:], kernel_size=15)
        
    
    new_soltab = solset.makeSoltab(soltype='phase', soltabName='XYoffset',
                             axesNames=['freq', 'ant', 'pol'], axesVals=[freq_per_sb, soltab.ant, ['XX','YY']],
                             vals=global_stat_offsets_smoothed,
                             weights=np.ones_like(global_stat_offsets_smoothed,dtype=np.float16))
    new_soltab.addHistory('CREATE (by PREFACTOR_XYOFFSET operation)')

    return 0


