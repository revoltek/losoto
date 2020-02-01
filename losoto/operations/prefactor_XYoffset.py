#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading PREFACTOR_XYOFFSET module.')

def _run_parser(soltab, parser, step):
    chanWidth = parser.getstr( step, 'chanWidth')

    parser.checkSpelling( step, soltab, ['chanWidth'])
    return run(soltab, chanWidth)


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


def run( soltab, chanWidth ):

    import numpy as np
    import scipy.signal as sg

    logging.info("Running XYoffset on: "+soltab.name)
    solset = soltab.getSolset()

    solType = soltab.getType()
    if solType != 'phase':
       logging.error("Soltab type of "+soltab.name+" is: "+solType+" but should be phase.")
       return 1

    for vals, coord, selection in soltab.getValuesIter(returnAxes=['pol','ant','freq', 'time'], weight=False):
        vals = reorderAxes( vals, soltab.getAxesNames(), ['time', 'ant', 'freq', 'pol'] )
        pass

    phases_tmp = np.copy(vals) # axes are [time, ant, freq, pol]
    freqs = np.copy(soltab.freq)
    npol = len(soltab.pol)
    refstationID=2

    stationsnames = [ stat for stat in soltab.ant]

    subbandHz = 195.3125e3
    if type(chanWidth) is str:
        letters = [1 for s in chanWidth[::-1] if s.isalpha()]
        indx = len(chanWidth) - sum(letters)
        unit = chanWidth[indx:]
        if unit.strip().lower() == 'hz':
            conversion = 1.0
        elif unit.strip().lower() == 'khz':
            conversion = 1e3
        elif unit.strip().lower() == 'mhz':
            conversion = 1e6
        else:
            logging.error("The unit on chanWidth was not understood.")
            raise ValueError("The unit on chanWidth was not understood.")
        chanWidthHz = float(chanWidth[:indx]) * conversion
    else:
        chanWidthHz = chanWidth
    offsetHz = subbandHz / 2.0 - 0.5 * chanWidthHz
    freqmin = np.min(soltab.freq[:]) + offsetHz # central frequency of first subband
    freqmax = np.max(soltab.freq[:]) + offsetHz # central frequency of last subband
    freqs_new  = np.arange(freqmin, freqmax+100e3, subbandHz)

    # this gets the subband number to any given frequency in HBA-low
    subbands = np.unique(np.round(freqs/195.3125e3-512.))
    nsubbands = len(subbands)
    nchan = len(freqs)/nsubbands
    if nsubbands*nchan != len(freqs):
        print("find_cal_global_phaseoffset_losoto.py: irregular number of ch/SB detected! Bailing out!")
        print("  nchan %d, nSB: %d, nfreq: %d" % (nchan, nsubbands, len(freqs)))
        sys.exit(1)
    tmpfreqs = freqs.reshape([nsubbands,nchan])
    freq_per_sb = np.mean(tmpfreqs,axis=1)
    nstations = len(stationsnames)
    refphases = phases_tmp[:, refstationID, :, :]

    for istat in range(nstations):
        phases_00 = phases_tmp[:, istat, :, 0] - refphases[:, :, 0]
        phases_11 = phases_tmp[:, istat, :, 1] - refphases[:, :, 1]
        phases_diff = normalize(phases_00 - phases_11)
        tmp_phases_diff = np.median(phases_diff, axis=0)  # take median over time axis
        med_phases_diff = np.median(tmp_phases_diff.reshape([nsubbands, nchan]), axis=1)  # take median over each subband
        if istat == 0:
            global_stat_offsets = med_phases_diff
        else:
            global_stat_offsets = np.vstack( (global_stat_offsets, med_phases_diff) )
    global_stat_offsets_smoothed = np.zeros([nsubbands, nstations, npol])
    global_stat_offsets_smoothed_interp = np.zeros([len(freqs_new), nstations, npol])
    for istat in range(nstations):
        global_stat_offsets_smoothed[:, istat, -1] = sg.medfilt(global_stat_offsets[istat, :], kernel_size=15) # smooth over frequency

        # Convert to real/imag, invert correction (so that offsets are removed when applied),
        # and interpolate to the output frequency grid
        real = np.interp(freqs_new, freq_per_sb, np.cos(-1. * global_stat_offsets_smoothed[:, istat, -1]))
        imag = np.interp(freqs_new, freq_per_sb, np.sin(-1. * global_stat_offsets_smoothed[:, istat, -1]))
        global_stat_offsets_smoothed_interp[:, istat, -1] = np.arctan2(imag, real)

    try:
        new_soltab = solset.getSoltab('XYoffset')
        new_soltab.delete()
    except:
        pass
    new_soltab = solset.makeSoltab(soltype='phase', soltabName='XYoffset',
                             axesNames=['freq', 'ant', 'pol'], axesVals=[freqs_new, soltab.ant, ['XX','YY']],
                             vals=global_stat_offsets_smoothed_interp,
                             weights=np.ones_like(global_stat_offsets_smoothed_interp, dtype=np.float16))
    new_soltab.addHistory('CREATE (by PREFACTOR_XYOFFSET operation)')

    return 0


