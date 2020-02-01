#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the direction-screen operation for LoSoTo

from losoto.lib_operations import *
from losoto.operations.stationscreen import _getxy, _radec2xy, _xy2radec, _makeWCS
from losoto.operations.stationscreen import _flag_outliers, _circ_chi2
from losoto._logging import logger as logging

logging.debug('Loading DIRECTIONSCREEN module.')


def _run_parser(soltab, parser, step):
    outSoltab = parser.getstr( step, "outSoltab", 'tecscreen' )
    height = parser.getfloat( step, "height", 200e3 )
    order = parser.getint( step, "Order", 5 )
    ncpu = parser.getint( '_global', "npcu", 0 )

    parser.checkSpelling( step, soltab, ['outSoltab', 'height', 'order'])
    return run(soltab, outSoltab, height, order, ncpu)


def _calculate_piercepoints(station_positions, source_positions, times, height=200e3):
    """
    Returns array of piercepoint locations and airmass values for a
    screen at the given height (in m)

    Parameters
    ----------
    station_positions : array
        array of station positions
    source_positions : array
        array of source positions
    times: array
        array of times
    height: float
        height of screen (m)

    Returns
    -------
    pp: array
        array of pierce points
    airmass: array
        array of airmass values
    midRA: float
        reference RA for WCS system (deg)
    midDec: float
        reference Dec for WCS system (deg)

    """
    import pyrap.measures
    import numpy as np
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar

    logging.info('Calculating screen pierce-point locations and airmass values...')
    N_sources = source_positions.shape[0]
    N_stations = station_positions.shape[0]
    N_piercepoints = N_stations * N_sources
    N_times = len(times)
    pp = np.zeros((N_times, N_piercepoints, 3))
    airmass = np.zeros((N_times, N_piercepoints))

    me = pyrap.measures.measures()
    position = me.position('ITRF', '%fm' % station_positions[0,0],
        '%fm' % station_positions[0,1], '%fm' % station_positions[0,2])
    me.doframe(position)

    pbar = progressbar.ProgressBar(maxval=N_times).start()
    ipbar = 0
    for k in range(N_times):
        epoch = me.epoch('UTC', '%fs' % times[k])
        me.doframe(epoch)
        pp_idx = 0
        for i in range(N_sources):
            ra = source_positions[i,0]
            dec = source_positions[i,1]
            d = me.direction('J2000', '%frad' % ra, '%frad' % dec)
            d1 = me.measure(d, 'ITRF')
            phi = d1['m0']['value']
            theta = d1['m1']['value']
            dx = np.cos(theta)*np.cos(phi)
            dy = np.cos(theta)*np.sin(phi)
            dz = np.sin(theta)
            direction = np.array((dx,dy,dz))
            for station_position in station_positions:
                pp[k, pp_idx, :], airmass[k, pp_idx] = _calc_piercepoint(
                    station_position, direction, height)
                pp_idx += 1
        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()
    midRA = 0.0
    midDec = 0.0

    return pp, airmass, midRA, midDec


def _calc_piercepoint(pos, direction, height):
    """
    Calculates pierce point locations and airmass values for given station
    position, source direction, and height (in m)

    Parameters
    ----------
    pos: array
        array of station positions
    direction: array
        array of source directions
    height: float
        height of screen (m)

    Returns
    -------
    pp: array
        array of pierce points
    airmass: array
        array of airmass values

    """
    import numpy as np

    pp = np.zeros(3)
    earth_ellipsoid_a = 6378137.0
    earth_ellipsoid_a2 = earth_ellipsoid_a * earth_ellipsoid_a
    earth_ellipsoid_b = 6356752.3142
    earth_ellipsoid_b2 = earth_ellipsoid_b * earth_ellipsoid_b
    earth_ellipsoid_e2 = (earth_ellipsoid_a2 - earth_ellipsoid_b2) / earth_ellipsoid_a2;

    ion_ellipsoid_a = earth_ellipsoid_a + height
    ion_ellipsoid_a2_inv = 1.0 / (ion_ellipsoid_a * ion_ellipsoid_a)
    ion_ellipsoid_b = earth_ellipsoid_b + height
    ion_ellipsoid_b2_inv = 1.0 / (ion_ellipsoid_b * ion_ellipsoid_b)

    x = pos[0] / ion_ellipsoid_a
    y = pos[1] / ion_ellipsoid_a
    z = pos[2] / ion_ellipsoid_b
    c = x*x + y*y + z*z - 1.0

    dx = direction[0] / ion_ellipsoid_a
    dy = direction[1] / ion_ellipsoid_a
    dz = direction[2] / ion_ellipsoid_b
    a = dx*dx + dy*dy + dz*dz
    b = x*dx + y*dy  + z*dz
    alpha = (-b + np.sqrt(b*b - a*c)) / a
    pp = pos[0] + alpha * direction
    normal_x = pp[0] * ion_ellipsoid_a2_inv
    normal_y = pp[1] * ion_ellipsoid_a2_inv
    normal_z = pp[2] * ion_ellipsoid_b2_inv
    norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z
    norm_normal = np.sqrt(norm_normal2)
    sin_lat2 = normal_z*normal_z / norm_normal2
    g = 1.0 - earth_ellipsoid_e2 * sin_lat2
    sqrt_g = np.sqrt(g)
    M = earth_ellipsoid_b2 / ( earth_ellipsoid_a * g * sqrt_g )
    N = earth_ellipsoid_a / sqrt_g

    local_ion_ellipsoid_e2 = (M - N) / ((M + height) * sin_lat2 - N - height);
    local_ion_ellipsoid_a = (N + height) * np.sqrt(1.0 - local_ion_ellipsoid_e2 * sin_lat2)
    local_ion_ellipsoid_b = local_ion_ellipsoid_a * np.sqrt(1.0 - local_ion_ellipsoid_e2)

    z_offset = ((1.0 - earth_ellipsoid_e2) * N + height - (1.0 -
        local_ion_ellipsoid_e2) * (N+height)) * np.sqrt(sin_lat2)

    x1 = pos[0] / local_ion_ellipsoid_a
    y1 = pos[1] / local_ion_ellipsoid_a
    z1 = (pos[2] - z_offset) / local_ion_ellipsoid_b
    c1 = x1*x1 + y1*y1 + z1*z1 - 1.0

    dx = direction[0] / local_ion_ellipsoid_a
    dy = direction[1] / local_ion_ellipsoid_a
    dz = direction[2] / local_ion_ellipsoid_b
    a = dx*dx + dy*dy + dz*dz
    b = x1*dx + y1*dy + z1*dz
    alpha = (-b + np.sqrt(b*b - a*c1)) / a

    pp = pos + alpha * direction

    normal_x = pp[0] * ion_ellipsoid_a2_inv
    normal_y = pp[1] * ion_ellipsoid_a2_inv
    normal_z = (pp[2] - z_offset) * ion_ellipsoid_b2_inv

    norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z
    norm_normal = np.sqrt(norm_normal2)

    airmass = norm_normal / (direction[0]*normal_x + direction[1]*normal_y +
        direction[2]*normal_z)

    return pp, airmass


def _fit_phase_screen(station_names, source_names, pp, airmass, rr, weights, times,
    height, order, r_0, beta, outQueue):
    """
    Fits a screen to given phase values using Karhunen-Lo`eve base vectors

    Parameters
    ----------
    station_names: array
        Array of station names
    source_names: array
        Array of source names
    pp: array
        Array of piercepoint locations
    airmass: array
        Array of airmass values (note: not currently used)
    rr: array
        Array of phase values to fit screen to
    weights: array
        Array of weights
    times: array
        Array of times
    height: float
        Height of screen (m)
    order: int
        Order of screen (i.e., number of KL base vectors to keep)
    r_0: float
        Scale size of phase fluctuations (m)
    beta: float
        Power-law index for phase structure function (5/3 => pure Kolmogorov
        turbulence)

    """
    import numpy as np
    from pylab import kron, concatenate, pinv, norm, newaxis, find, amin, svd, eye

    logging.info('Fitting screens...')

    # Initialize arrays
    N_stations = len(station_names)
    N_sources = len(source_names)
    N_times = len(times)
    N_piercepoints = N_sources * N_stations
    real_fit_white_all = np.zeros((N_times, N_sources, N_stations))
    imag_fit_white_all = np.zeros((N_times, N_sources, N_stations))
    phase_fit_white_all = np.zeros((N_times, N_sources, N_stations))
    real_residual_all = np.zeros((N_times, N_sources, N_stations))
    imag_residual_all = np.zeros((N_times, N_sources, N_stations))
    phase_residual_all = np.zeros((N_times, N_sources, N_stations))

    # Change phase to real/imag
    rr_real = np.cos(rr)
    rr_imag = np.sin(rr)

    for k in range(N_times):
        try:
            D = np.resize(pp[k, :, :], (N_piercepoints, N_piercepoints, 3))
            D = np.transpose(D, (1, 0, 2)) - D
            D2 = np.sum(D**2, axis=2)
            C = -(D2 / r_0**2)**(beta / 2.0) / 2.0
            pinvC = pinv(C, rcond=1e-3)
            U, S, V = svd(C)
            invU = pinv(np.dot(np.transpose(U[:, :order]), np.dot(weights[:, :, k], U[:, :order])), rcond=1e-3)

            # Calculate real screen
            rr1 = np.dot(np.transpose(U[:, :order]), np.dot(weights[:, :, k], rr_real[:, k]))
            real_fit = np.dot(pinvC, np.dot(U[:, :order], np.dot(invU, rr1)))
            real_fit_white_all[k, :, :] = real_fit.reshape((N_sources, N_stations))
            residual = rr_real - np.dot(C, real_fit)[:, newaxis]
            real_residual_all[k, :, :] = residual.reshape((N_sources, N_stations))

            # Calculate imag screen
            rr1 = np.dot(np.transpose(U[:, :order]), np.dot(weights[:, :, k], rr_imag[:, k]))
            imag_fit = np.dot(pinvC, np.dot(U[:, :order], np.dot(invU, rr1)))
            imag_fit_white_all[k, :, :] = imag_fit.reshape((N_sources, N_stations))
            residual = rr_imag - np.dot(C, imag_fit)[:, newaxis]
            imag_residual_all[k, :, :] = residual.reshape((N_sources, N_stations))

            # Calculate phase screen
            phase_fit = np.dot(pinvC, np.arctan2(np.dot(C, imag_fit), np.dot(C, real_fit)))
            phase_fit_white_all[k, :, :] = phase_fit.reshape((N_sources, N_stations))
            residual = rr - np.dot(C, phase_fit)[:, newaxis]
            phase_residual_all[k, :, :] = residual.reshape((N_sources, N_stations))
        except:
            # Set screen to zero if fit did not work
            logging.debug('Screen fit failed for timeslot {}'.format(k))
            real_fit_white_all[k, :, :] = np.zeros((N_sources, N_stations))
            real_residual_all[k, :, :] = np.ones((N_sources, N_stations))
            imag_fit_white_all[k, :, :] = np.zeros((N_sources, N_stations))
            imag_residual_all[k, :, :] = np.ones((N_sources, N_stations))
            phase_fit_white_all[k, :, :] = np.zeros((N_sources, N_stations))
            phase_residual_all[k, :, :] = np.ones((N_sources, N_stations))

    outQueue.put([real_fit_white_all, real_residual_all,
                  imag_fit_white_all, imag_residual_all,
                  phase_fit_white_all, phase_residual_all,
                  times])


def _fit_tec_screen(station_names, source_names, pp, airmass, rr, weights, times,
    height, order, r_0, beta, outQueue):
    """
    Fits a screen to given TEC values using Karhunen-Lo`eve base vectors

    Parameters
    ----------
    station_names: array
        Array of station names
    source_names: array
        Array of source names
    pp: array
        Array of piercepoint locations
    airmass: array
        Array of airmass values (note: not currently used)
    rr: array
        Array of TEC values to fit screen to
    weights: array
        Array of weights
    times: array
        Array of times
    height: float
        Height of screen (m)
    order: int
        Order of screen (i.e., number of KL base vectors to keep)
    r_0: float
        Scale size of phase fluctuations (m)
    beta: float
        Power-law index for phase structure function (5/3 => pure Kolmogorov
        turbulence)

    """
    import numpy as np
    from pylab import kron, concatenate, pinv, norm, newaxis, find, amin, svd, eye

    logging.info('Fitting screens...')

    # Initialize arrays
    N_stations = len(station_names)
    N_sources = len(source_names)
    N_times = len(times)
    N_piercepoints = N_sources * N_stations
    tec_fit_white_all = np.zeros((N_times, N_sources, N_stations))
    tec_residual_all = np.zeros((N_times, N_sources, N_stations))

    for k in range(N_times):
        D = np.resize(pp[k, :, :], (N_piercepoints, N_piercepoints, 3))
        D = np.transpose(D, (1, 0, 2)) - D
        D2 = np.sum(D**2, axis=2)
        C = -(D2 / r_0**2)**(beta / 2.0) / 2.0
        pinvC = pinv(C, rcond=1e-3)
        U, S, V = svd(C)
        invU = pinv(np.dot(np.transpose(U[:, :order]), np.dot(weights[:, :, k], U[:, :order])), rcond=1e-3)

        # Calculate screen
        rr1 = np.dot(np.transpose(U[:, :order]), np.dot(weights[:, :, k], rr[:, k]))
        tec_fit = np.dot(pinvC, np.dot(U[:, :order], np.dot(invU, rr1)))
        tec_fit_white_all[k, :, :] = tec_fit.reshape((N_sources, N_stations))
        residual = rr - np.dot(C, tec_fit)[:, newaxis]
        tec_residual_all[k, :, :] = residual.reshape((N_sources, N_stations))

    outQueue.put([tec_fit_white_all, tec_residual_all, times])


def run(soltab, outSoltab='tecscreen', height=200.0e3, order=12,
    beta=5.0/3.0, ncpu=0):
    """
    Fits a screen to TEC + scalaraphase values.

    The results of the fit are stored in the soltab parent solset in
    "outSoltab" and the residual phases (actual-screen) are stored in
    "outsoltabresid". These values are the screen phase values per station per
    pierce point per solution interval. The pierce point locations are stored in
    an auxiliary array in the output soltabs.

    Screens can be plotted with the PLOTSCREEN operation.

    Parameters
    ----------
    soltab: solution table
        Soltab containing phase solutions
    outSoltab: str, optional
        Name of output soltab
    height : float, optional
        Height in m of screen
    order : int, optional
        Order of screen (i.e., number of KL base vectors to keep).
    beta: float, optional
        Power-law index for phase structure function (5/3 => pure Kolmogorov
        turbulence)
    ncpu: int, optional
        Number of CPUs to use. If 0, all are used
    niter: int, optional
        Number of iterations to do when determining weights
    nsigma: float, optional
        Number of sigma above which directions are flagged

    """
    import numpy as np
    from numpy import newaxis
    import re
    import os

    # Get screen type
    screen_type = soltab.getType()
    if screen_type not in ['phase', 'tec']:
        logging.error('Screens can only be fit to soltabs of type "phase" or "tec".')
        return 1
    logging.info('Using solution table {0} to calculate {1} screens'.format(soltab.name, screen_type))

    # Load phases
    axis_names = soltab.getAxesNames()
    r = np.array(soltab.val)
    weights = soltab.weight[:]
    if 'freq' in soltab.getAxesNames():
        freqs = soltab.freq[:]
        if len(freqs) > 1:
            logging.error('Screens can only be fit at a single frequency')
            return 1
        freq = freqs[0]

        # remove degenerate freq axis
        freq_ind = soltab.getAxesNames().index('freq')
        r = np.squeeze(r, axis=freq_ind)
        weights = np.squeeze(weights, axis=freq_ind)
        axis_names.pop(freq_ind)

    # fix for missing dir axis
    if not 'dir' in soltab.getAxesNames():
        r = np.array([r])
        weights = np.array([weights])
        dir_ind = len(axis_names)
        source_names = ['POINTING']
    else:
        dir_ind = axis_names.index('dir')
        source_names = soltab.dir[:]

    time_ind = axis_names.index('time')
    ant_ind = axis_names.index('ant')
    r = r.transpose([dir_ind, ant_ind, time_ind])
    weights = weights.transpose([dir_ind, ant_ind, time_ind])
    times = np.array(soltab.time)

    # Collect station and source names and positions and times, making sure
    # that they are ordered correctly.
    solset = soltab.getSolset()
    source_dict = solset.getSou()
    source_positions = []
    for source in source_names:
        source_positions.append(source_dict[source])
    station_names = soltab.ant[:]
    station_dict = solset.getAnt()
    station_positions = []
    for station in station_names:
        station_positions.append(station_dict[station])
    N_sources = len(source_names)
    N_times = len(times)
    N_stations = len(station_names)

    logging.info('Using height = {0} m and order = {1}'.format(height, order))
    if height < 100e3:
        logging.warning("Height is less than 100e3 m. This is likely too low.")

    # Initialize various arrays
    N_piercepoints = N_sources * N_stations
    if screen_type == 'phase':
        real_screen = np.zeros((N_sources, N_stations, N_times))
        real_residual =np.zeros((N_sources, N_stations, N_times))
        imag_screen = np.zeros((N_sources, N_stations, N_times))
        imag_residual = np.zeros((N_sources, N_stations, N_times))
    screen = np.zeros((N_sources, N_stations, N_times))
    residual = np.zeros((N_sources, N_stations, N_times))
    val_amp = 1.0
    r_0 = 100.0 # shouldn't matter what we choose
    rr = np.reshape(r, [N_piercepoints, N_times])

    # Find pierce points and airmass values for given screen height
    pp, airmass, midRA, midDec = _calculate_piercepoints(np.array(station_positions),
        np.array(source_positions), np.array(times), height)

    # Fit the screens
    station_weights = np.reshape(weights, [N_piercepoints, N_times])
    if screen_type == 'phase':
        mpm = multiprocManager(ncpu, _fit_phase_screen)
        for tindx, t in enumerate(times):
            w = np.diag(station_weights[:, tindx])[:, :, newaxis]
            mpm.put([station_names, source_names, pp[tindx, newaxis, :, :],
                airmass[tindx, newaxis, :], rr[:, tindx, newaxis], w,
                [t], height, order, r_0, beta])
        mpm.wait()
        for (real_scr, real_res, imag_scr, imag_res, phase_scr, phase_res, t) in mpm.get():
            i = times.tolist().index(t[0])
            real_screen[:, :, i] = real_scr[0, :, :]
            real_residual[:, :, i] = real_res[0, :, :]
            imag_screen[:, :, i] = imag_scr[0, :, :]
            imag_residual[:, :, i] = imag_res[0, :, :]
            screen[:, :, i] = phase_scr[0, :, :]
            residual[:, :, i] = phase_res[0, :, :]
    elif screen_type == 'tec':
        mpm = multiprocManager(ncpu, _fit_tec_screen)
        for tindx, t in enumerate(times):
            w = np.diag(station_weights[:, tindx])[:, :, newaxis]
            mpm.put([station_names, source_names, pp[tindx, newaxis, :, :],
                airmass[tindx, newaxis, :], rr[:, tindx, newaxis], w,
                [t], height, order, r_0, beta])
        mpm.wait()
        for (scr, res, t) in mpm.get():
            i = times.tolist().index(t[0])
            screen[:, :, i] = scr[0, :, :]
            residual[:, :, i] = res[0, :, :]
    weights = np.reshape(station_weights, (N_sources, N_stations, N_times))

    # Write the results to the output solset
    dirs_out = source_names
    times_out = times
    ants_out = station_names

    # Store screen values
    weights = weights.transpose([0, 2, 1]) # order is now [dir, time, ant]
    vals = screen.transpose([0, 2, 1])
    screen_st = solset.makeSoltab('{}screen'.format(screen_type), outSoltab,
        axesNames=['dir', 'time', 'ant'], axesVals=[dirs_out, times_out,
        ants_out], vals=vals, weights=weights)
    vals = residual.transpose([0, 2, 1])
    resscreen_st = solset.makeSoltab('{}screenresid'.format(screen_type), outSoltab+'resid',
        axesNames=['dir', 'time', 'ant'], axesVals=[dirs_out, times_out,
        ants_out], vals=vals, weights=weights)

    # Store beta, r_0, height, and order as attributes of the screen soltabs
    screen_st.obj._v_attrs['beta'] = beta
    screen_st.obj._v_attrs['r_0'] = r_0
    screen_st.obj._v_attrs['height'] = height
    screen_st.obj._v_attrs['order'] = order
    if 'freq' in soltab.getAxesNames():
        screen_st.obj._v_attrs['freq'] = freq

    # Store piercepoint table. Note that it does not conform to the axis
    # shapes, so we cannot use makeSoltab()
    solset.obj._v_file.create_array('/'+solset.name+'/'+screen_st.obj._v_name,
        'piercepoint', obj=pp)

    screen_st.addHistory('CREATE (by DIRECTIONSCREEN operation)')
    resscreen_st.addHistory('CREATE (by DIRECTIONSCREEN operation)')

    return 0
