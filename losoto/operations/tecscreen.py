#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the TEC-screen operation for LoSoTo

import logging
from losoto.operations_lib import *

logging.debug('Loading TECSCREEN module.')


def calculate_piercepoints(station_positions, source_positions, times, height = 200e3):
    """
    Returns array of piercepoint locations and airmass values for a
    screen at the given height (in m)

    Keyword arguments:
    station_positions -- array of station positions
    source_positions -- array of source positions
    times -- array of times
    height -- height of screen (m)
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

    me = pyrap.measures.measures()
    position = me.position('ITRF', '%fm' % station_positions[0,0],
        '%fm' % station_positions[0,1], '%fm' % station_positions[0,2])
    me.doframe(position)

    pp = np.zeros((N_times, N_piercepoints,3))
    airmass = np.zeros((N_times, N_piercepoints))

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

    return pp, airmass


def calc_piercepoint(pos, direction, height):
    """
    Calculates pierce point locations and airmass values for given station
    position, source direction, and height (in m)

    Keyword arguments:
    pos -- array of station positions
    direction -- array of source directions
    height -- height of screen (m)
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
    b = x1*dx + y1*dy  + z1*dz
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


def fit_screen_to_tec(station_names, source_names, pp, airmass, rr, times,
    height, order, r_0, beta):
    """
    Fits a screen to given TEC values using Karhunen-Lo`eve base vectors

    Keyword arguments:
    station_names -- array of station names
    source_names -- array of source names
    pp -- array of piercepoint locations
    airmass -- array of airmass values
    rr -- array of TEC solutions
    times -- array of times
    height -- height of screen (m)
    order -- order of screen (i.e., number of KL base vectors to keep)
    r_0 -- scale size of phase fluctuations (m)
    beta -- power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    """
    import numpy as np
    from pylab import kron, concatenate, pinv, norm, newaxis, find, amin, svd, eye
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar

    logging.info('Fitting screens to TEC values...')
    N_stations = len(station_names)
    N_sources = len(source_names)
    N_times = len(times)

    tec_fit_all = np.zeros((N_times, N_sources, N_stations))
    residual_all = np.zeros((N_times, N_sources, N_stations))

    A = concatenate([kron(eye(N_sources), np.ones((N_stations, 1))),
        kron(np.ones((N_sources, 1)), eye(N_stations))], axis=1)

    N_piercepoints = N_sources * N_stations
    P = eye(N_piercepoints) - np.dot(np.dot(A, pinv(np.dot(A.T, A))), A.T)

    pbar = progressbar.ProgressBar(maxval=N_times).start()
    ipbar = 0
    for k in range(N_times):
        try:
            D = np.resize(pp[k, :, :], (N_piercepoints, N_piercepoints, 3))
            D = np.transpose(D, (1, 0, 2)) - D
            D2 = np.sum(D**2, axis=2)
            C = -(D2 / r_0**2)**(beta / 2.0) / 2.0
            P1 = eye(N_piercepoints) - np.ones((N_piercepoints, N_piercepoints)) / N_piercepoints
            C1 = np.dot(np.dot(P1, C), P1)
            U, S, V = svd(C1)

            B = np.dot(P, np.dot(np.diag(airmass[k, :]), U[:, :order]))
            pinvB = pinv(B, rcond=1e-3)

            rr1 = np.dot(P, rr[:, k])
            tec_fit = np.dot(U[:, :order], np.dot(pinvB, rr1))
            tec_fit_all[k, :, :] = tec_fit.reshape((N_sources, N_stations))

            residual = rr1 - np.dot(P, tec_fit)
            residual_all[k, :, :] = residual.reshape((N_sources, N_stations))
        except:
            # Set screen to zero if fit did not work
            logging.debug('Tecscreen fit failed for timeslot {0}'.format(k))
            tec_fit_all[k, :, :] = np.zeros((N_sources, N_stations))
            residual_all[k, :, :] = np.ones((N_sources, N_stations))

        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()

    return tec_fit_all, residual_all


def run( step, parset, H ):
    """
    Fits a screen to TEC values derived by the TECFIT operation.

    The TEC values are read from the specified tec soltab.

    The results of the fit are stored in the specified tecscreen solution table.
    These values are the screen TEC values per station per pierce point per
    solution interval. The pierce point locations are stored in an auxiliary
    array in the output solution table.

    TEC screens can be plotted with the PLOT operation by setting PlotType =
    TECScreen.

    The H5parm_exporter.py tool can be used to export the screen to a parmdb
    that BBS and the AWimager can use. Note, however, that the output screens
    are not normalized properly (any normalization was lost due to the use of
    source-to-source phase gradients in the TECFIT operation). Therefore, a
    direction-independent calibration must be done after exporting the screens
    to a parmdb file, with the following settings in the BBS solve step:

        Model.Ionosphere.Enable = T
        Model.Ionosphere.Type = EXPION
    """
    import numpy as np
    import re
    from losoto.h5parm import solFetcher, solWriter
    # Switch to the Agg backend to prevent problems with pylab imports when
    # DISPLAY env. variable is not set
    import os

    soltabs = getParSoltabs( step, parset, H )
    outSoltabs = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "OutSoltab"]), [] )
    height = np.array(parset.getDoubleVector('.'.join(["LoSoTo.Steps", step, "Height"]), [200e3] ))
    order = int(parset.getString('.'.join(["LoSoTo.Steps", step, "Order"]), '5' ))

    # Load TEC values from TECFIT operation
    indx = 0
    for soltab in openSoltabs(H, soltabs):
        if 'tec' not in soltab._v_title:
            logging.warning('No TECFIT solution tables found for solution table '
                '{0}'.format(soltabs[indx]))
            continue
            indx += 1
        solset = soltabs[indx].split('/')[0]
        logging.info('Using input solution table: {0}'.format(soltabs[indx]))
        logging.info('Using output solution table: {0}'.format(outSoltabs[indx]))

        # Collect station and source names and positions and times, making sure
        # that they are ordered correctly.
        t = solFetcher(soltab)
        r, axis_vals = t.getValues()
        try:
            source_names = axis_vals['dir']
        except:
            source_names = ['pointing']
        source_dict = H.getSou(solset)
        source_positions = []
        for source in source_names:
            source_positions.append(source_dict[source])
        station_names = axis_vals['ant']
        station_dict = H.getAnt(solset)
        station_positions = []
        for station in station_names:
            station_positions.append(station_dict[station])
        times = axis_vals['time']

        # Get sizes
        N_sources = len(source_names)
        N_times = len(times)
        N_stations = len(station_names)
        N_piercepoints = N_sources * N_stations
        if N_sources > 1:
            if len(r.shape) > 3: r = np.squeeze(r) # remove degenerate freq axis added by ndppp
            # I want: dir ant time
            rr = np.reshape(r.transpose([1, 0, 2]), [N_piercepoints, N_times])
        else:
            assert t.getAxesNames()[0] == 'time'
            assert t.getAxesNames()[1] == 'ant'
            rr = np.reshape(r.transpose([1,0]), [N_piercepoints, N_times])

        heights = list(set(np.linspace(height[0], height[-1], 5)))
        heights.sort()
        if len(heights) > 1:
            logging.info('Trying range of heights: {0} m'.format(heights))
        for i, height in enumerate(heights):
            # Find pierce points and airmass values for given screen height
            logging.info('Using height = {0} m and order = {1}'.format(height, order))
            if height < 100e3:
                logging.warning("Height is less than 100e3 m.")
            pp, airmass = calculate_piercepoints(np.array(station_positions),
                np.array(source_positions), np.array(times), height)

            # Fit a TEC screen
            r_0 = 10e3
            beta = 5.0 / 3.0
            tec_screen, residual = fit_screen_to_tec(station_names, source_names,
                pp, airmass, rr, times, height, order, r_0, beta)
            total_resid = np.sum(np.abs(residual))
            if i > 0:
                if total_resid < best_resid:
                    tec_screen_best = tec_screen
                    pp_best = pp
                    height_best = height
                    best_resid = total_resid
            else:
                tec_screen_best = tec_screen
                pp_best = pp
                height_best = height
                best_resid = total_resid
            if len(heights) > 1:
                logging.info('Total residual for fit: {0}\n'.format(total_resid))

        # Use screen with lowest total residual
        if len(heights) > 1:
            tec_screen = tec_screen_best
            pp = pp_best
            height = height_best
            logging.info('Using height (with lowest total residual) of {0} m'.format(height))

        # Write the results to the output solset
        dirs_out = source_names
        times_out = times
        ants_out = station_names

        # Make output tecscreen table
        outSolset = outSoltabs[indx].split('/')[0]
        outSoltab = outSoltabs[indx].split('/')[1]
        if not outSolset in H.getSolsets().keys():
            solsetTEC = H.makeSolset(outSolset)
            dirs_pos = source_positions
            sourceTable = solsetTEC._f_get_child('source')
            sourceTable.append(zip(*(dirs_out, dirs_pos)))
            ants_pos = station_positions
            antennaTable = solsetTEC._f_get_child('antenna')
            antennaTable.append(zip(*(ants_out, ants_pos)))

        # Store tecscreen values. The residual values are stored in the weights
        # table. Flagged values of the screen have weights set to 0.0.
        vals = tec_screen.transpose([1, 0, 2])
        weights = residual.transpose([1, 0, 2])
        tec_screen_st = H.makeSoltab(outSolset, 'tecscreen', outSoltab,
            axesNames=['dir', 'time', 'ant'], axesVals=[dirs_out, times_out,
            ants_out], vals=vals, weights=weights)

        # Store beta, r_0, height, and order as attributes of the tecscreen
        # soltab
        tec_screen_st._v_attrs['beta'] = beta
        tec_screen_st._v_attrs['r_0'] = r_0
        tec_screen_st._v_attrs['height'] = height
        tec_screen_st._v_attrs['order'] = order

        # Make output piercepoint table
        tec_screen_solset = tec_screen_st._v_parent._v_name
        H.H.create_carray('/'+tec_screen_solset+'/'+tec_screen_st._v_name,
            'piercepoint', obj=pp)

        # Add histories
        sw = solWriter(tec_screen_st)
        sw.addHistory('CREATE (by TECSCREEN operation)')
        indx += 1

    return 0
