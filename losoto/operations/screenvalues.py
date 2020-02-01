#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the operation to get values from screens

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading SCREENVALUES module.')


def _run_parser(soltab, parser, step):
    inSoltab1 = parser.getstr( step, "inSoltab1" )
    inSoltab2 = parser.getstr( step, "inSoltab2", None )
    outSoltab = parser.getstr( step, "outSoltab" )
    sourceDict = parser.getstr( step, "sourceDict" )
    ncpu = parser.getint( step, "ncpu", 0 )

    parser.checkSpelling( step, soltab, ['inSoltab1', 'sourceDict', 'outSoltab', 'inSoltab2'])
    return run(inSoltab1, sourceDict, outSoltab, inSoltab2, ncpu)

def _calculate_tecsp(screen1, screen2, pp, directions, k, sindx, beta_val,
    r_0, freq1, freq2, midRA, midDec, outQueue):
    """
    Calculates TEC and scalar phase from screens at two frequencies

    screen1: array
        Array of screen values at the piercepoints for freq1
    screen2: array
        Array of screen values at the piercepoints for freq2
    pp: array
        Array of piercepoint locations
    directions: array
        Array of RA, Dec in degrees for which phases are desired
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
    import numpy as np
    from losoto.operations.stationscreen import _radec2xy

    # Convert direction (RA, Dec) to (x, y)
    x, y = _radec2xy(directions[0], directions[1], midRA, midDec)

    # Calculate phases at freq1 and freq2
    phase1 = np.zeros(len(x))
    phase2 = np.zeros(len(x))
    for i, (xi, yi) in enumerate(zip(x, y)):
        p = np.array([xi, yi, 0.0])
        d2 = np.sum(np.square(pp - p), axis=1)
        c = -(d2 / ( r_0**2 ))**(beta_val / 2.0) / 2.0
        phase1[i] = np.dot(c, screen1)
        phase2[i] = np.dot(c, screen2)

    # Calculate TEC and scalarphase from phase1 and phase2
    tec = normalize_phase(phase1 - phase2) / (1.0/freq1 - 1.0/freq2) / (-8.4479745e9)
    scalarphase = normalize_phase((freq1*phase1 - freq2*phase2) / (freq1 - freq2))

    outQueue.put([k, tec, scalarphase])


def _calculate_val(screen, pp, directions, k, sindx, beta_val, r_0, midRA,
    midDec, outQueue):
    """
    Calculates values from screen

    screen: array
        Array of screen values at the piercepoints
    pp: array
        Array of piercepoint locations
    directions: array
        Array of RA, Dec in degrees for which phases are desired
    r_0: float
        Scale size of phase fluctuations
    beta_val: float
        Power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees

    """
    import numpy as np
    from losoto.operations.stationscreen import _radec2xy

    # Convert direction (RA, Dec) to (x, y)
    x, y = _radec2xy(directions[0], directions[1], midRA, midDec)

    # Calculate values
    N_dirs = len(x)
    N_freqs = screen.shape[1]
    values = np.zeros((N_dirs, N_freqs))
    for i, (xi, yi) in enumerate(zip(x, y)):
        p = np.array([xi, yi, 0.0])
        d2 = np.sum(np.square(pp - p), axis=1)
        c = -(d2 / ( r_0**2 ))**(beta_val / 2.0) / 2.0
        for findx in range(N_freqs):
            values[i, findx] = np.dot(c, screen[:, findx])

    outQueue.put([k, values])


def _screens_to_tecsp(pp, screen1, screen2, directions, station_positions,
    beta_val, r_0, freq1, freq2, midRA=0.0, midDec=0.0, ncpu=0):
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
    height: float
        Height of screen (m)
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
    import numpy as np
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar

    ra = [r for r, d in directions]
    dec = [d for r, d in directions]
    directions = np.array([ra, dec]) * 180.0 / np.pi

    N_sources = len(ra)
    N_times = screen1.shape[1]
    N_freqs = 1
    N_stations = screen1.shape[2]
    tec = np.zeros((N_times, N_stations, N_sources, N_freqs))
    scalarphase = np.zeros((N_times, N_stations, N_sources, N_freqs))

    logging.info('Calculating phases...')
    N_total = N_stations
    pbar = progressbar.ProgressBar(maxval=N_total).start()
    ipbar = 0
    for sindx in range(N_stations):
        mpm = multiprocManager(ncpu, _calculate_tecsp)
        for tindx in range(N_times):
            mpm.put([screen1[:, tindx, sindx], screen2[:, tindx, sindx],
                pp, directions, tindx, sindx, beta_val, r_0,
                freq1, freq2, midRA, midDec])
        mpm.wait()
        for (tindx, t, s) in mpm.get():
            tec[tindx, sindx, :, 0] = t
            scalarphase[tindx, sindx, :, 0] = s
        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()

    return (tec, scalarphase)


def _screen_to_val(pp, screen, directions, station_positions, beta_val, r_0,
    midRA=0.0, midDec=0.0, ncpu=0):
    """
    Caculates phases from screens

    Parameters
    ----------
    pp: array
        Array of piercepoint locations
    screen: array
        Array of screen values at the piercepoints with axes [dir_ind, time_ind,
        freq_ind, ant_ind, pol_ind]
    directions: list of arrays
        List of (RA, Dec) arrays in radians for which phases are desired
    height: float
        Height of screen (m)
    r_0: float
        Scale size of phase fluctuations
    beta_val: float
        Power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees
    ncpu: int
        Number of CPUs to use

    """
    import numpy as np
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar

    ra = [r for r, d in directions]
    dec = [d for r, d in directions]
    directions = np.array([ra, dec]) * 180.0 / np.pi

    N_sources = len(ra)
    N_times = screen.shape[1]
    N_freqs = screen.shape[2]
    N_stations = screen.shape[3]
    if len(screen.shape) == 4:
        N_pols = 1
        values = np.zeros((N_times, N_stations, N_sources, N_freqs))
    else:
        N_pols = screen.shape[4]
        values = np.zeros((N_times, N_stations, N_sources, N_freqs, N_pols))

    logging.info('Calculating values...')
    N_total = N_pols * N_stations
    pbar = progressbar.ProgressBar(maxval=N_total).start()
    ipbar = 0
    for sindx in range(N_stations):
        for pindx in range(N_pols):
            mpm = multiprocManager(ncpu, _calculate_val)
            for tindx in range(N_times):
                if len(screen.shape) == 4:
                    inscreen = screen[:, tindx, :, sindx]
                else:
                    inscreen = screen[:, tindx, :, sindx, pindx]
                mpm.put([inscreen, pp, directions, tindx, sindx,
                    beta_val, r_0, midRA, midDec])
            mpm.wait()
            for (tindx, val) in mpm.get():
                if len(screen.shape) == 4:
                    values[tindx, sindx, :, :] = val
                else:
                    values[tindx, sindx, :, :, pindx] = val
            pbar.update(ipbar)
            ipbar += 1
    pbar.finish()

    return values


def run(soltab1, source_dict, outsoltab, soltab2=None, ncpu=0):
    """
    Calculate phases from screens

    Parameters
    ----------
    soltab1 : solution table
        Soltab containing the screen at frequency 1
    source_dict : dict
        Dictionary of source positions for which phases are desired
    outsoltab : str
        Name of output solution table. If TEC+scalarphase solutions are desired,
        two solution tables are made: 'tecoutsoltab' and 'scalarphaseoutsoltab'
    soltab2 : solution table, optional
        Soltab containing the screens at frequency 2 (needed for TEC+scalar
        phase)
    ncpu : int, optional
        Number of CPUs to use. If 0, all are used

    """
    import numpy as np

    if soltab2 is not None:
        screen_type = 'tecsp'
        logging.info('Using input solution tables {0} and {1} to calculate '
            'TEC+scalarphase values'.format(soltab1.name, soltab2.name))
    else:
        screen_type = soltab1.getType().split('screen')[0]
        if screen_type not in ['phase', 'amplitude', 'tec']:
            logging.error('Values can only be derived for screens of type "phasescreen", "tecscreen", or "amplitudescreen".')
            return 1
        logging.info('Using input solution table {0} to calculate '
            '{1} values'.format(soltab1.name, screen_type))

    # Get values and reorder axes
    solset = soltab1.getSolset()
    screen1 = np.array(soltab1.val)
    times = np.array(soltab1.time)
    if screen_type == 'tecsp':
        screen2 = np.array(soltab2.val)
    axis_names = soltab1.getAxesNames()
    freq_ind = axis_names.index('freq')
    dir_ind = axis_names.index('dir')
    time_ind = axis_names.index('time')
    ant_ind = axis_names.index('ant')
    if 'pol' in axis_names:
        pol_ind = axis_names.index('pol')
        screen1 = screen1.transpose([dir_ind, time_ind, freq_ind, ant_ind, pol_ind])
    else:
        screen1 = screen1.transpose([dir_ind, time_ind, freq_ind, ant_ind])
        if screen_type == 'tecsp':
            screen2 = screen2.transpose([dir_ind, time_ind, freq_ind, ant_ind])
    if screen_type == 'tecsp':
        # For TEC+scalarphase solutions, the two screens must
        # both be at a single (identical) frequency, so remove the freq axes
        screen1 = np.squeeze(screen1, axis=2)
        screen2 = np.squeeze(screen2, axis=2)


    # Collect info
    source_names = list(source_dict.keys())
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
    beta_val = soltab1.obj._v_attrs['beta']
    r_0 = soltab1.obj._v_attrs['r_0']
    if screen_type == 'tecsp':
        freq1 = soltab1.freq[0]
        freq2 = soltab2.freq[0]
        if freq1 == freq2:
            logging.error('Screens have the same frequency: cannot calculate '
                'TEC and scalarphase values!'.format(f))
            return 1
    pp = soltab1.obj.piercepoint[:]
    midRA = soltab1.obj._v_attrs['midra']
    midDec = soltab1.obj._v_attrs['middec']

    if screen_type == 'tecsp':
        tecs, scalarphases = _screens_to_tecsp(pp, screen1, screen2,
            source_positions, np.array(station_positions), beta_val, r_0,
            freq1, freq2, midRA=midRA, midDec=midDec, ncpu=ncpu)
    else:
        values = _screen_to_val(pp, screen1, source_positions,
            np.array(station_positions), beta_val, r_0, midRA=midRA,
            midDec=midDec, ncpu=ncpu)

    # Update solset source table if new sources are present
    needs_update = False
    for k, v in source_dict.items():
        if k not in solset.getSou():
            needs_update = True
            break
    if needs_update:
        source_dict_orig = solset.getSou()
        solset.obj.source.remove()
        descriptor = np.dtype([('name', np.str_, 128),('dir', np.float32, 2)])
        soltab = solset.obj._v_file.create_table(solset.obj, 'source', descriptor, title = 'Source names and directions', expectedrows = 25)
        sourceTable = solset.obj._f_get_child('source')
        names = []
        positions = []
        for k, v in source_dict_orig.items():
            # Add sources from original dict
            names.append(k)
            if type(v) is list:
                positions.append(v)
            else:
                positions.append(v.tolist())
        for k, v in source_dict.items():
            # Add any new sources
            if k not in names:
                names.append(k)
                if type(v) is list:
                    positions.append(v)
                else:
                    positions.append(v.tolist())
        sourceTable.append(list(zip(*(names, positions))))

    # Write the results to the output soltab(s)
    dirs_out = source_names
    times_out = times
    ants_out = station_names
    if screen_type == 'tecsp':
        freqs_out = [(freq1 + freq2) / 2.0]
        weights = np.ones(tecs.shape)
        st1 = solset.makeSoltab('tec', 'tec'+outsoltab, axesNames=['time',
            'ant', 'dir', 'freq'], axesVals=[times_out, ants_out, dirs_out,
            freqs_out], vals=tecs, weights=weights)
        st2 = solset.makeSoltab('scalarphase', 'scalarphase'+outsoltab,
            axesNames=['time', 'ant', 'dir', 'freq'], axesVals=[times_out,
            ants_out, dirs_out, freqs_out], vals=scalarphases,
            weights=weights)
    else:
        freqs_out = soltab1.freq[:]
        weights = np.ones(values.shape)
        if 'pol' in axis_names:
            pols_out = soltab1.pol[:]
            st = solset.makeSoltab(screen_type, outsoltab, axesNames=['time',
                'ant', 'dir', 'freq', 'pol'], axesVals=[times_out, ants_out,
                dirs_out, freqs_out, pols_out], vals=values, weights=weights)
        else:
            st = solset.makeSoltab(screen_type, outsoltab, axesNames=['time',
                'ant', 'dir', 'freq'], axesVals=[times_out, ants_out,
                dirs_out, freqs_out], vals=values, weights=weights)

    return 0
