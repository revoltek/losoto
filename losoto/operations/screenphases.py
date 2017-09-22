#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the screen-to-phases operation for LoSoTo


import logging
from losoto.operations_lib import *

logging.debug('Loading SCREENPHASES module.')


def run_parser(soltab, parser, step):
    inSoltab1 = parser.getarraystr( step, "inSoltab1" )
    inSoltab2 = parser.getarraystr( step, "inSoltab2" )
    outSoltab = parser.getarraystr( step, "outSoltab" )

    return run(inSoltab1, inSoltab2, sourceDict, frequencies, outSoltab, ncpu)


def calculate_phases(screen1, screen2, pp, directions, frequencies,
    N_piercepoints, k, sindx, beta_val, r_0, freq1, freq2, midRA, midDec,
    outQueue):
    """
    Calculates phases from screens

    screen1: array
        Array of screen values at the piercepoints for freq1
    screen2: array
        Array of screen values at the piercepoints for freq2
    pp: array
        Array of piercepoint locations
    directions: array
        Array of RA, Dec in degrees for which phases are desired
    frequencies: array
        Array of frequencies in Hz for which phases are desired
    order: int
        Order of screen (e.g., number of KL base vectors to keep)
    r_0: float
        Scale size of phase fluctuations
    beta_val: float
        Power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    freq1: float
        Frequency of screen1 in Hz
    freq2: float
        Frequency of screen2 in Hz
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees

    """
    from numpy import kron, concatenate, newaxis
    from numpy.linalg import pinv, norm
    import numpy as np
    from losoto.operations.tecscreen import calc_piercepoint
    from losoto.operations.phasescreen import radec2xy

    # Convert direction (RA, Dec) to (x, y)
    x, y = radec2xy(directions[0], directions[1], midRA, midDec)

    # Calculate phases at freq1 and freq2
    f1 = screen1.reshape(N_piercepoints)
    f2 = screen2.reshape(N_piercepoints)
    phase1 = np.zeros(len(x))
    phase2 = np.zeros(len(x))
    for i, (xi, yi) in enumerate(zip(x, y)):
        p = np.array([xi, yi, 0.0])
        d2 = np.sum(np.square(pp - p), axis=1)
        c = -(d2 / ( r_0**2 ))**(beta_val / 2.0) / 2.0
        phase1[i] = np.dot(c, f1)
        phase2[i] = np.dot(c, f2)

    # Interpolate phases to desired frequencies:
    # TEC = (ph1 - ph2) / (1/freq1 - 1/freq2) / -8.4479745e9
    # scph = (freq1*ph1 - freq2*ph2) / (freq1 - freq2)
    phase = np.zeros((len(x), len(frequencies)))
    for i, freq in enumerate(frequencies):
        tec = (phase1 - phase2) / (1.0/freq1 - 1.0/freq2) / (-8.4479745e9)
        scphase = (freq1*phase1 - freq2*phase2) / (freq1 - freq2)
        phase[:, i] = -8.4479745e9 * tec/freq + scphase # order is [dir, freq]

    outQueue.put([k, phase])


def screen_to_phases(pp, screen1, screen2, directions, frequencies,
    station_positions, times, order, beta_val, r_0,
    freq1, freq2, midRA=0.0, midDec=0.0, ncpu=0):
    """
    Caculates phases from screens

    Parameters
    ----------
    pp: array
        Array of piercepoint locations
    screen1: array
        Array of screen values at the piercepoints for freq1
    screen2: array
        Array of screen values at the piercepoints for freq2
    directions: list of arrays
        List of (RA, Dec) arrays in radians for which phases are desired
    frequencies: array
        Array of frequencies in Hz for which phases are desired
    times: array
        Array of times
    height: float
        Height of screen (m)
    order: int
        Order of screen (e.g., number of KL base vectors to keep)
    r_0: float
        Scale size of phase fluctuations
    beta_val: float
        Power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    freq1: float
        Frequency of screen1 in Hz
    freq2: float
        Frequency of screen2 in Hz
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees
    ncpu: int
        Number of CPUs to use

    """
    from numpy import kron, concatenate, newaxis
    from numpy.linalg import pinv, norm
    import numpy as np
    import os
    from losoto.operations.tecscreen import calc_piercepoint
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar

    # input check
    if ncpu == 0:
        import multiprocessing
        ncpu = multiprocessing.cpu_count()

    ra = [r for r, d in directions]
    dec = [d for r, d in directions]
    directions = np.array([ra, dec]) * 180.0 / np.pi

    N_sources = len(ra)
    N_times = len(times)
    N_freqs = len(frequencies)
    N_stations = station_positions.shape[0]
    phase = np.zeros((N_sources, N_stations, N_freqs, N_times))
    for sindx in range(station_positions.shape[0]):
        N_piercepoints = N_sources

        logging.info('Calculating phases...')
        mpm = multiprocManager(ncpu, calculate_phases)
        for k in range(N_times):
            mpm.put([screen1[:, k, sindx], screen2[:, k, sindx], pp[k, :, :],
                directions, frequencies, N_piercepoints, k, sindx, beta_val,
                r_0, freq1, freq2, midRA, midDec])
        mpm.wait()
        for (k, ph) in mpm.get():
            phase[:, sindx, :, k] = ph # order is [dir, ant, freq, time]

    return phase


def run(soltab1, soltab2, source_dict, frequencies, outsoltab, ncpu=0):
    """
    Plot screens (one plot is made per time and per station)

    Parameters
    ----------
    soltab1 : solution table
        Soltab containing the screen at frequency 1
    soltab2: solution table
        Soltab containing the screen at frequency 2
    source_dict: dict
        Dictionary of source positions for which phases are desired
    frequencies: array
        Array of frequencies in Hz for which phases are desired
    outsoltab: str
        Name of output solution table
    ncpu: int, optional
        Number of CPUs to use. If 0, all are used

    """
    import os
    import numpy as np

    logging.info('Using input solution tables: {0} and {1}'.format(soltab1.name,
        soltab2.name))

    # Get values from soltabs
    solset = soltab1.getSolset()
    screen1 = np.array(soltab1.val)
    screen2 = np.array(soltab2.val)
    times = np.array(soltab1.time)

    # Collect info
    source_names = source_dict.keys()
    source_positions = []
    for source in source_names:
        source_positions.append(source_dict[source])
    station_names = soltab1.ant
    station_dict = solset.getAnt()
    station_positions = []
    for station in station_names:
        station_positions.append(station_dict[station])
    height = soltab1.obj._v_attrs['height']
    if height != 0.0:
        logging.error('Screen must have height = 0 to use this operation')
        return 1
    order = soltab1.obj._v_attrs['order']
    beta_val = soltab1.obj._v_attrs['beta']
    r_0 = soltab1.obj._v_attrs['r_0']
    freq1 = soltab1.obj._v_attrs['freq']
    freq2 = soltab2.obj._v_attrs['freq']
    pp = soltab1.obj.piercepoint
    midRA = soltab1.obj._v_attrs['midra']
    midDec = soltab1.obj._v_attrs['middec']

    phases = screen_to_phases(pp, screen1, screen2, source_positions,
        frequencies, np.array(station_positions), times, order, beta_val,
        r_0, freq1, freq2, midRA=midRA, midDec=midDec, ncpu=ncpu)

    # Write the results to the output soltab
    dirs_out = source_names
    times_out = times
    ants_out = station_names
    freqs_out = frequencies
    weights = np.ones(phases.shape)
    st = solset.makeSoltab(outsoltab, axesNames=['dir', 'ant', 'freq', 'time'],
        axesVals=[dirs_out, ants_out, freqs_out, times_out], vals=phases,
        weights=weights)

    return 0
