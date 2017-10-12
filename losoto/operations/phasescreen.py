#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the phase-screen operation for LoSoTo


import logging
from losoto.operations_lib import *

logging.debug('Loading PHASESCREEN module.')


def run_parser(soltab, parser, step):
    outSoltab = parser.getarraystr( step, "outSoltab" )
    height = parser.getfloat( step, "height", 200e3 )
    order = parser.getint( step, "Order", 5 )
    return run(soltab, outSoltab, height, order)


def calculate_piercepoints(station_positions, source_positions, times, height=200e3):
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
        height of screen (m). Set to 0.0 to fit to RA, Dec instead of plane

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

    if height == 0.0:
        xyz = np.zeros((N_sources, 3))
        ra_deg = source_positions.T[0] * 180.0 / np.pi
        dec_deg = source_positions.T[1] * 180.0 / np.pi
        xy, midRA, midDec = getxy(ra_deg, dec_deg)
        xyz[:, 0] = xy[0]
        xyz[:, 1] = xy[1]
        for k in range(N_times):
            pp_idx = 0
            for i in range(N_sources):
                for station_position in station_positions:
                    pp[k, pp_idx, :] = xyz[i]
                    pp_idx += 1
    else:
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
                    pp[k, pp_idx, :], airmass[k, pp_idx] = calc_piercepoint(
                        station_position, direction, height)
                    pp_idx += 1
            pbar.update(ipbar)
            ipbar += 1
        pbar.finish()
        midRA = 0.0
        midDec = 0.0

    return pp, airmass, midRA, midDec


def calc_piercepoint(pos, direction, height):
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


def getxy(RA, Dec, midRA=None, midDec=None):
    """
    Returns array of projected x and y values.

    Parameters
    ----------
    RA, Dec : list
        Lists of RA and Dec in degrees
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees

    Returns
    -------
    x, y : numpy array, numpy array, float, float
        arrays of x and y values

    """
    import numpy as np

    if midRA is None or midDec is None:
        x, y  = radec2xy(RA, Dec)

        # Refine x and y using midpoint
        if len(x) > 1:
            xmid = min(x) + (max(x) - min(x)) / 2.0
            ymid = min(y) + (max(y) - min(y)) / 2.0
            xind = np.argsort(x)
            yind = np.argsort(y)
            try:
                midxind = np.where(np.array(x)[xind] > xmid)[0][0]
                midyind = np.where(np.array(y)[yind] > ymid)[0][0]
                midRA = RA[xind[midxind]]
                midDec = Dec[yind[midyind]]
                x, y  = radec2xy(RA, Dec, midRA, midDec)
            except IndexError:
                midRA = RA[0]
                midDec = Dec[0]
        else:
            midRA = RA[0]
            midDec = Dec[0]

    x, y  = radec2xy(RA, Dec, refRA=midRA, refDec=midDec)

    return np.array([x, y]), midRA, midDec


def radec2xy(RA, Dec, refRA=None, refDec=None):
    """
    Returns x, y for input ra, dec.

    Note that the reference RA and Dec must be the same in calls to both
    radec2xy() and xy2radec() if matched pairs of (x, y) <=> (RA, Dec) are
    desired.

    Parameters
    ----------
    RA : list
        List of RA values in degrees
    Dec : list
        List of Dec values in degrees
    refRA : float, optional
        Reference RA in degrees.
    refDec : float, optional
        Reference Dec in degrees

    Returns
    -------
    x, y : list, list
        Lists of x and y pixel values corresponding to the input RA and Dec
        values

    """
    import numpy as np

    x = []
    y = []
    if refRA is None:
        refRA = RA[0]
    if refDec is None:
        refDec = Dec[0]

    # Make wcs object to handle transformation from ra and dec to pixel coords.
    w = makeWCS(refRA, refDec)

    for ra_deg, dec_deg in zip(RA, Dec):
        ra_dec = np.array([[ra_deg, dec_deg]])
        x.append(w.wcs_world2pix(ra_dec, 0)[0][0])
        y.append(w.wcs_world2pix(ra_dec, 0)[0][1])

    return x, y


def xy2radec(x, y, refRA=0.0, refDec=0.0):
    """
    Returns x, y for input ra, dec.

    Note that the reference RA and Dec must be the same in calls to both
    radec2xy() and xy2radec() if matched pairs of (x, y) <=> (RA, Dec) are
    desired.

    Parameters
    ----------
    x : list
        List of x values in pixels
    y : list
        List of y values in pixels
    refRA : float, optional
        Reference RA in degrees
    refDec : float, optional
        Reference Dec in degrees

    Returns
    -------
    RA, Dec : list, list
        Lists of RA and Dec values corresponding to the input x and y pixel
        values

    """
    import numpy as np

    RA = []
    Dec = []

    # Make wcs object to handle transformation from ra and dec to pixel coords.
    w = makeWCS(refRA, refDec)

    for xp, yp in zip(x, y):
        x_y = np.array([[xp, yp]])
        RA.append(w.wcs_pix2world(x_y, 0)[0][0])
        Dec.append(w.wcs_pix2world(x_y, 0)[0][1])

    return RA, Dec


def makeWCS(refRA, refDec):
    """
    Makes simple WCS object.

    Parameters
    ----------
    refRA : float
        Reference RA in degrees
    refDec : float
        Reference Dec in degrees

    Returns
    -------
    w : astropy.wcs.WCS object
        A simple TAN-projection WCS object for specified reference position

    """
    from astropy.wcs import WCS
    import numpy as np

    w = WCS(naxis=2)
    w.wcs.crpix = [1000, 1000]
    w.wcs.cdelt = np.array([-0.0005, 0.0005])
    w.wcs.crval = [refRA, refDec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.set_pv([(2, 1, 45.0)])

    return w


def normalize(phase):
    """
    Normalize phase (in rad) to the range [-pi, pi].
    """
    import numpy as np

    # Convert to range [-2*pi, 2*pi].
    out = np.fmod(phase, 2.0 * np.pi)

    # Convert to range [-pi, pi]
    out[out < -np.pi] += 2.0 * np.pi
    out[out > np.pi] -= 2.0 * np.pi

    return out


def estimate_weights_median_window(srcindx, phases, Nmedian, Nstddev, outQueue):
    """
    Set weights using a median-filter method

    """
    import numpy as np
    from losoto.operations.phasescreen import normalize

    # Change phase to real/imag
    real = np.cos(phases)
    imag = np.sin(phases)

    # Smooth with short window and subtract
    N = Nmedian
    idx = np.arange(N) + np.arange(len(phases)-N+1)[:, None]
    med = np.median(real[idx], axis=1)
    med_real = np.zeros(phases.shape)
    med_real[0:(N-1)/2] = med[0]
    med_real[(N-1)/2:-(N-1)/2] = med
    med_real[-(N-1)/2:] = med[-1]
    med = np.median(imag[idx], axis=1)
    med_imag = np.zeros(phases.shape)
    med_imag[0:(N-1)/2] = med[0]
    med_imag[(N-1)/2:-(N-1)/2] = med
    med_imag[-(N-1)/2:] = med[-1]
    diff_real = real - med_real
    diff_imag = imag - med_imag
    c = diff_real + diff_imag*1j

    # Calculate standard deviation in larger window
    N = Nstddev
    idx = np.arange(N) + np.arange(len(phases)-N+1)[:, None]
    mad = np.std(c[idx], axis=1)
    mad_c = np.zeros(phases.shape)
    mad_c[0:(N-1)/2] = mad[0]
    mad_c[(N-1)/2:-(N-1)/2] = mad
    mad_c[-(N-1)/2:] = mad[-1]
    nzeros = np.where(c == 0.0)[0].shape[0]
    corr_factor = 2.0 * float(phases.shape[0]) / float(nzeros) # compensate for reduction in scatter due to smoothing
    weights = 1.0 / np.square(corr_factor*mad_c)

    outQueue.put([srcindx, weights])


def flag_outliers(srcindx, rr, weights, phase_residual, nsigma, N, screen_stddev, outQueue):
    """
    Flags outliers

    Parameters
    ----------
    srcindx: int
        Index of direction
    rr: array
        Array of actual phases
    weights: array
        Array of weights
    phase_residual: array
        Array of residual values from phase screen fitting (rad)
    nsigma: float
        Number of sigma above with outliers are clipped (= weight set to zero)
    N: int
        Number of time samples to use for sliding window
    screen_stddev: float
        Screen standard dev over all directions

    Returns
    -------
    srcindx: int
        Index of direction
    weights: array
        array of weights


    """
    import numpy as np

    # Compare smoothed residuals to stddev of station and screen
    stddev = np.sqrt(1.0/weights)
    phase_residual = normalize(phase_residual)
    outlier_ind = np.where(np.logical_or((phase_residual > nsigma*screen_stddev),
        (phase_residual > nsigma*stddev)))
    weights[outlier_ind] = 0.0

    outQueue.put([srcindx, weights])


def fit_phase_screen(station_names, source_names, pp, airmass, rr, weights, times,
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


def run(soltab, outsoltab='phasescreen', height=0.0, order=12,
    beta=5.0/3.0, ncpu=0, niter=3, nsigma=3.0):
    """
    Fits a screen to TEC + scalaraphase values.

    The results of the fit are stored in the soltab parent solset in
    "outsoltab" and the residual phases (actual-screen) are stored in
    "outsoltabresid". These values are the screen phase values per station per
    pierce point per solution interval. The pierce point locations are stored in
    an auxiliary array in the output soltabs.

    Screens can be plotted with the PLOTSCREEN operation.

    Parameters
    ----------
    soltab: solution table
        Soltab containing phase solutions
    outsoltab: str, optional
        Name of output soltab
    height : float, optional
        Height in m of screen. If 0, the RA and Dec values (projected on the
        image plane) are used instead for the screen pierce points
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

    # input check
    if ncpu == 0:
        import multiprocessing
        ncpu = multiprocessing.cpu_count()

    # Load phases
    logging.info('Using solution table: {}'.format(soltab.name))
    r = np.array(soltab.val)
    weights = soltab.weight[:]
    if len(r.shape) > 3:
        # remove degenerate freq axis
        freq_ind = soltab.getAxesNames().index('freq')
        r = np.squeeze(r, axis=freq_ind)
        weights = np.squeeze(weights, axis=freq_ind)
    r = r.transpose([2, 1, 0]) # order is now [dir, ant, time]
    weights = weights.transpose([2, 1, 0])
    times = np.array(soltab.time)
    freqs = soltab.freq[:]
    if len(freqs) > 1:
        logging.error('Screens can only be fit at a single frequency')
        return 1
    freq = freqs[0]

    # Collect station and source names and positions and times, making sure
    # that they are ordered correctly.
    solset = soltab.getSolset()
    source_names = soltab.dir[:]
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
    if height == 0.0:
        # Fit station screens
        logging.info("Using RA and Dec for pierce points")

        # Initialize various arrays
        real_screen = np.zeros((N_sources, N_stations, N_times))
        real_residual =np.zeros((N_sources, N_stations, N_times))
        imag_screen = np.zeros((N_sources, N_stations, N_times))
        imag_residual = np.zeros((N_sources, N_stations, N_times))
        phase_screen = np.zeros((N_sources, N_stations, N_times))
        phase_residual = np.zeros((N_sources, N_stations, N_times))

        # Fit screens
        val_amp = 1.0
        r_0 = 100 # shouldn't matter what we choose
        for s, stat in enumerate(station_names):

            N_piercepoints = N_sources
            rr = np.reshape(r[:, s, :], [N_piercepoints, N_times])

            # Find pierce points and airmass values for given screen height
            pp, airmass, midRA, midDec = calculate_piercepoints(np.array([station_positions[s]]),
                np.array(source_positions), np.array(times), height)

            # Iterate:
            # 2. flag nsigma outliers
            # 3. refit with new weights
            # 4. repeat for niter
            station_weights = weights[:, s, :]
            init_station_weights = weights[:, s, :].copy() # preserve initial weights
            for iterindx in range(niter):
                if iterindx > 0:
                    # Flag outliers
                    nonflagged = np.where(init_station_weights > 0.0)
                    real_diff = np.cos(rr[nonflagged]) - np.cos(rr[nonflagged]-phase_residual[:, s, :][nonflagged])
                    imag_diff = np.sin(rr[nonflagged]) - np.sin(rr[nonflagged]-phase_residual[:, s, :][nonflagged])
                    stddev_real = np.sqrt(np.average(real_diff**2, weights=init_station_weights[nonflagged], axis=0))
                    stddev_imag = np.sqrt(np.average(imag_diff**2, weights=init_station_weights[nonflagged], axis=0))
                    total_stddev = np.arctan2(stddev_real, stddev_imag)
                    mpm = multiprocManager(ncpu, flag_outliers)
                    for srcindx in range(N_piercepoints):
                        mpm.put([srcindx, rr[srcindx, :], init_station_weights[srcindx, :],
                            phase_residual[srcindx, s, :], nsigma, 3, total_stddev])
                    mpm.wait()
                    for (srcindx, w) in mpm.get():
                        station_weights[srcindx, :] = w

                # Fit the screens
                mpm = multiprocManager(ncpu, fit_phase_screen)
                for tindx, t in enumerate(times):
                    w = np.diag(station_weights[:, tindx])[:, :, newaxis]
                    mpm.put([[stat], source_names, pp[tindx, newaxis, :, :],
                        airmass[tindx, newaxis, :], rr[:, tindx, newaxis], w,
                        [t], height, order, r_0, beta])
                mpm.wait()
                for (real_scr, real_res, imag_scr, imag_res, phase_scr, phase_res, t) in mpm.get():
                    i = times.tolist().index(t[0])
                    real_screen[:, s, i] = real_scr[0, :, 0]
                    real_residual[:, s, i] = real_res[0, :, 0]
                    imag_screen[:, s, i] = imag_scr[0, :, 0]
                    imag_residual[:, s, i] = imag_res[0, :, 0]
                    phase_screen[:, s, i] = phase_scr[0, :, 0]
                    phase_residual[:, s, i] = phase_res[0, :, 0]
            weights[:, s, :] = station_weights

    else:
        # Fit global screens
        if height < 100e3:
            logging.warning("Height is less than 100e3 m. This is likely too low.")

        # Initialize various arrays
        N_piercepoints = N_sources * N_stations
        real_screen = np.zeros((N_sources, N_stations, N_times))
        real_residual =np.zeros((N_sources, N_stations, N_times))
        imag_screen = np.zeros((N_sources, N_stations, N_times))
        imag_residual = np.zeros((N_sources, N_stations, N_times))
        phase_screen = np.zeros((N_sources, N_stations, N_times))
        phase_residual = np.zeros((N_sources, N_stations, N_times))

        # Fit screens
        val_amp = 1.0
        r_0 = 100.0 # shouldn't matter what we choose
        rr = np.reshape(r, [N_piercepoints, N_times])

        # Find pierce points and airmass values for given screen height
        pp, airmass, midRA, midDec = calculate_piercepoints(np.array(station_positions),
            np.array(source_positions), np.array(times), height)

        # Iterate:
        # 2. flag nsigma outliers
        # 3. refit with new weights
        # 4. repeat for niter
        station_weights = np.reshape(weights, [N_piercepoints, N_times])
        init_station_weights = station_weights.copy() # preserve initial weights
        for iterindx in range(niter):
            if iterindx > 0:
                # Flag outliers
                phase_residual_r = np.reshape(phase_residual, [N_piercepoints, N_times])
                nonflagged = np.where(init_station_weights > 0.0)
                real_diff = np.cos(rr[nonflagged]) - np.cos(rr[nonflagged]-phase_residual_r[nonflagged])
                imag_diff = np.sin(rr[nonflagged]) - np.sin(rr[nonflagged]-phase_residual_r[nonflagged])
                stddev_real = np.std(real_diff, axis=0)
                stddev_imag = np.std(imag_diff, axis=0)
                total_stddev = np.arctan2(stddev_real, stddev_imag)
                mpm = multiprocManager(ncpu, flag_outliers)
                for srcindx in range(N_piercepoints):
                    mpm.put([srcindx, rr[srcindx, :], init_station_weights[srcindx, :],
                        phase_residual_r[srcindx, :], nsigma, 3, total_stddev])
                mpm.wait()
                for (srcindx, w) in mpm.get():
                    station_weights[srcindx, :] = w

            # Fit the screens
            mpm = multiprocManager(ncpu, fit_phase_screen)
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
                phase_screen[:, :, i] = phase_scr[0, :, :]
                phase_residual[:, :, i] = phase_res[0, :, :]
        weights = np.reshape(station_weights, (N_sources, N_stations, N_times))

    # Write the results to the output solset
    dirs_out = source_names
    times_out = times
    ants_out = station_names

    # Store tecscreen values
    weights = weights.transpose([0, 2, 1]) # order is now [dir, time, ant]
    vals = phase_screen.transpose([0, 2, 1])
    screen_st = solset.makeSoltab('phasescreen', outsoltab,
        axesNames=['dir', 'time', 'ant'], axesVals=[dirs_out, times_out,
        ants_out], vals=vals, weights=weights)
    vals = phase_residual.transpose([0, 2, 1])
    resscreen_st = solset.makeSoltab('phasescreenresid', outsoltab+'resid',
        axesNames=['dir', 'time', 'ant'], axesVals=[dirs_out, times_out,
        ants_out], vals=vals, weights=weights)

    # Store beta, r_0, height, and order as attributes of the screen soltabs
    screen_st.obj._v_attrs['beta'] = beta
    screen_st.obj._v_attrs['r_0'] = r_0
    screen_st.obj._v_attrs['height'] = height
    screen_st.obj._v_attrs['freq'] = freq
    screen_st.obj._v_attrs['order'] = order
    if height == 0.0:
        screen_st.obj._v_attrs['midra'] = midRA
        screen_st.obj._v_attrs['middec'] = midDec

    # Store piercepoint table. Note that it does not conform to the axis
    # shapes, so we cannot use makeSoltab()
    solset.obj._v_file.create_array('/'+solset.name+'/'+screen_st.obj._v_name,
        'piercepoint', obj=pp)

    screen_st.addHistory('CREATE (by PHASESCREEN operation)')
    resscreen_st.addHistory('CREATE (by PHASESCREEN operation)')

    return 0
