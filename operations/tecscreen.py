#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the TEC-screen operation for LoSoTo

import logging
from operations_lib import *
import pyrap.measures
import numpy
import os
from pylab import *


logging.debug('Loading TECSCREEN module.')


def calculate_piercepoints(station_positions, source_positions, times, height = 200e3):
    """
    Returns array of piercepoint locations and airmass values for a
    screen at the given height
    """
    import progressbar

    logging.info('Calculating pierce points...')
    N_sources = source_positions.shape[0]
    N_stations = station_positions.shape[0]
    N_piercepoints = N_stations * N_sources
    N_times = len(times)

    me = pyrap.measures.measures()
    position = me.position('ITRF', '%fm' % station_positions[0,0], '%fm' % station_positions[0,1], '%fm' % station_positions[0,2])
    me.doframe(position)


    pp = numpy.zeros((N_times, N_piercepoints,3))
    airmass = numpy.zeros((N_times, N_piercepoints))

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
            dx = numpy.cos(theta)*numpy.cos(phi)
            dy = numpy.cos(theta)*numpy.sin(phi)
            dz = numpy.sin(theta)
            direction = numpy.array((dx,dy,dz))
            for station_position in station_positions:
                pp[k, pp_idx, :], airmass[k, pp_idx] = calc_piercepoint(station_position, direction, height)
                pp_idx += 1
        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()

    return pp, airmass


def calc_piercepoint(pos, direction, height):
   pp = zeros(3)
   earth_ellipsoid_a = 6378137.0;
   earth_ellipsoid_a2 = earth_ellipsoid_a*earth_ellipsoid_a;
   earth_ellipsoid_b = 6356752.3142;
   earth_ellipsoid_b2 = earth_ellipsoid_b*earth_ellipsoid_b;
   earth_ellipsoid_e2 = (earth_ellipsoid_a2 - earth_ellipsoid_b2) / earth_ellipsoid_a2;

   ion_ellipsoid_a = earth_ellipsoid_a + height;
   ion_ellipsoid_a2_inv = 1.0 / (ion_ellipsoid_a * ion_ellipsoid_a);
   ion_ellipsoid_b = earth_ellipsoid_b + height;
   ion_ellipsoid_b2_inv = 1.0 / (ion_ellipsoid_b * ion_ellipsoid_b);

   x = pos[0]/ion_ellipsoid_a;
   y = pos[1]/ion_ellipsoid_a;
   z = pos[2]/ion_ellipsoid_b;
   c = x*x + y*y + z*z - 1.0;

   dx = direction[0] / ion_ellipsoid_a
   dy = direction[1] / ion_ellipsoid_a
   dz = direction[2] / ion_ellipsoid_b
   a = dx*dx + dy*dy + dz*dz
   b = x*dx + y*dy  + z*dz
   alpha = (-b + sqrt(b*b - a*c))/a
   pp = pos[0] + alpha*direction
   normal_x = pp[0] * ion_ellipsoid_a2_inv
   normal_y = pp[1] * ion_ellipsoid_a2_inv
   normal_z = pp[2] * ion_ellipsoid_b2_inv
   norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z
   norm_normal = sqrt(norm_normal2)
   sin_lat2 = normal_z*normal_z / norm_normal2
   g = 1.0 - earth_ellipsoid_e2*sin_lat2
   sqrt_g = sqrt(g)
   M = earth_ellipsoid_b2 / ( earth_ellipsoid_a * g * sqrt_g )
   N = earth_ellipsoid_a / sqrt_g

   local_ion_ellipsoid_e2 = (M-N) / ((M+height)*sin_lat2 - N - height);
   local_ion_ellipsoid_a = (N+height) * sqrt(1.0 - local_ion_ellipsoid_e2*sin_lat2)
   local_ion_ellipsoid_b = local_ion_ellipsoid_a*sqrt(1.0 - local_ion_ellipsoid_e2)

   z_offset = ((1.0-earth_ellipsoid_e2)*N + height - (1.0-local_ion_ellipsoid_e2)*(N+height)) * sqrt(sin_lat2)

   x1 = pos[0]/local_ion_ellipsoid_a
   y1 = pos[1]/local_ion_ellipsoid_a
   z1 = (pos[2]-z_offset)/local_ion_ellipsoid_b
   c1 = x1*x1 + y1*y1 + z1*z1 - 1.0

   dx = direction[0] / local_ion_ellipsoid_a
   dy = direction[1] / local_ion_ellipsoid_a
   dz = direction[2] / local_ion_ellipsoid_b
   a = dx*dx + dy*dy + dz*dz
   b = x1*dx + y1*dy  + z1*dz
   alpha = (-b + sqrt(b*b - a*c1))/a

   pp = pos + alpha*direction

   normal_x = pp[0] * ion_ellipsoid_a2_inv
   normal_y = pp[1] * ion_ellipsoid_a2_inv
   normal_z = (pp[2] - z_offset) * ion_ellipsoid_b2_inv

   norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z
   norm_normal = sqrt(norm_normal2)

   airmass = norm_normal / (direction[0]*normal_x + direction[1]*normal_y + direction[2]*normal_z)

   return pp, airmass


def fit_screen_to_tec(station_names, source_names, pp, airmass, rr, times,
    height, order, r_0, beta):
    """
    Fits a screen to given TEC values
    """
    import progressbar

    logging.info('Fitting screens to TEC values...')
    N_stations = len(station_names)
    N_sources = len(source_names)
    N_times = len(times)

    tec_fit_all = zeros((N_times, N_sources, N_stations))
    tec_fit_white_all = zeros((N_times, N_sources, N_stations))

    A = concatenate([kron(eye(N_sources), ones((N_stations,1))),
        kron(ones((N_sources,1)), eye(N_stations))], axis=1)

    N_piercepoints = N_sources*N_stations
    P = eye(N_piercepoints) - dot(dot(A, pinv(dot(A.T, A))), A.T)

    pbar = progressbar.ProgressBar(maxval=N_times).start()
    ipbar = 0
    for k in range(N_times):
        D = resize( pp[k,:,:], ( N_piercepoints, N_piercepoints, 3 ) )
        D = transpose( D, ( 1, 0, 2 ) ) - D
        D2 = sum( D**2, axis=2 )
        C = -(D2 / ( r_0**2 ) )**( beta / 2.0 )/2.0
        P1 = eye(N_piercepoints) - ones((N_piercepoints, N_piercepoints)) / N_piercepoints
        C1 = dot(dot(P1, C ), P1)
        U,S,V = svd(C1)

        B = dot(P, dot(diag(airmass[k,:]), U[:,:order]))
        pinvB = pinv(B, rcond=1e-3)

        rr1 = dot(P, rr[:,k])

        tec_fit = dot(U[:, :order], dot(pinvB, rr1))
        tec_fit_all[k, :, :] = tec_fit.reshape((N_sources, N_stations))
        tec_fit_white = dot(pinv(C), tec_fit)
        tec_fit_white_all[k, :, :] = tec_fit_white.reshape((N_sources, N_stations))
        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()

    return tec_fit_all, tec_fit_white_all


def run( step, parset, H ):
    """
    Fits a screen to TEC values derived by the TECFIT operation.

    The TEC values are read from the specified tec-type soltab.

    The results of the fit are stored in the specified tecfitwhite- and
    piercepoint-type soltabs.

    TEC screens can be plotted with the PLOT operation by setting PlotType =
    TECScreen.

    Note that the output screens are not normalized (any normalization was lost
    due to the use of source-to-source phase gradients in the TECFIT operation).
    Therefore, a direction-independent calibration must be done after exporting
    the screens to a parmdb file, with the following settings in the BBS solve
    step:
        Model.Ionosphere.Enable = T
        Model.Ionosphere.Type = EXPION
    """
    import numpy as np
    import re
    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    outSoltabsTEC = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "OutSoltabTEC"]), [] )
    outSoltabsPP = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "OutSoltabPP"]), [] )
    height = np.float(parset.getString('.'.join(["LoSoTo.Steps", step, "Height"]), '200' ))
    order = int(parset.getString('.'.join(["LoSoTo.Steps", step, "Order"]), '15' ))
    logging.info('Using height = {0} m and order = {1}'.format(height, order))

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
        logging.info('Using output solution tables: {0} and {1}'.format(
            outSoltabsTEC[indx], outSoltabsPP[indx]))
        station_dict = H.getAnt(solset)
        station_names = station_dict.keys()
        station_positions = station_dict.values()
        source_dict = H.getSou(solset)
        source_names = source_dict.keys()
        source_positions = source_dict.values()

        t = solFetcher(soltab)
        r, axis_vals = t.getValues()
        times = axis_vals['time']
        N_sources = len(source_names)
        N_times = len(times)
        N_stations = len(station_names)
        N_piercepoints = N_sources * N_stations
        rr = reshape(r.transpose([0,2,1]), [N_piercepoints, N_times])

        # Find pierce points and airmass values for given screen height
        pp, airmass = calculate_piercepoints(np.array(station_positions),
            np.array(source_positions), np.array(times), height)

        # Fit a TEC screen
        r_0 = 10e3
        beta = 5.0/3.0
        tec_fit, tec_fit_white = fit_screen_to_tec(station_names, source_names,
            pp, airmass, rr, times, height, order, r_0, beta)

        # Write the results to the output solset
        dirs_out = []
        dirs_pos = []
        for s in range(N_sources):
            dirs_out.append(source_names[s])

        times_out = times

        ants_out = []
        ants_pos = []
        for s in range(N_stations):
            ants_out.append(station_names[s])

        # Make output tecfitwhite table
        outSolsetTEC = outSoltabsTEC[indx].split('/')[0]
        outSoltabTEC = outSoltabsTEC[indx].split('/')[1]
        if not outSolsetTEC in H.getSolsets().keys():
            solsetTEC = H.makeSolset(outSolsetTEC)
            dirs_out = []
            dirs_pos = []
            for s in range(N_sources):
                dirs_out.append(source_names[s])
                dirs_pos.append(source_positions[s])
            sourceTable = solsetTEC._f_get_child('source')
            sourceTable.append(zip(*(dirs_out, dirs_pos)))
            ants_out = []
            ants_pos = []
            for s in range(N_stations):
                ants_out.append(station_names[s])
                ants_pos.append(station_positions[s])
            antennaTable = solsetTEC._f_get_child('antenna')
            antennaTable.append(zip(*(ants_out, ants_pos)))

        # Store tecfitwhite values. The tec_fit values are stored in the weights
        # table to simplify things
        vals = tec_fit_white.transpose([1, 0, 2])
        weights = tec_fit.transpose([1, 0, 2])
        tec_fit_st = H.makeSoltab(outSolsetTEC, 'tecfitwhite', outSoltabTEC,
            axesNames=['dir', 'time', 'ant'], axesVals=[dirs_out, times_out,
            ants_out], vals=vals, weights=weights)

        # Store beta, r_0, height, and order as attributes of the tecfitwhite
        # soltab
        tec_fit_st._v_attrs['beta'] = beta
        tec_fit_st._v_attrs['r_0'] = r_0
        tec_fit_st._v_attrs['height'] = height
        tec_fit_st._v_attrs['order'] = order

        # Make output piercepoint table
        outSolsetPP = outSoltabsPP[indx].split('/')[0]
        outSoltabPP = outSoltabsPP[indx].split('/')[1]
        if not outSolsetPP in H.getSolsets().keys():
            solsetPP = H.makeSolset(outSolsetPP)
            sourceTable = solsetPP._f_get_child('source')
            sourceTable.append(zip(*(dirs_out, dirs_pos)))
            antennaTable = solsetPP._f_get_child('antenna')
            antennaTable.append(zip(*(ants_out, ants_pos)))

        # Store piercepoint positions
        pp_indx = range(N_piercepoints)
        pp_pos_indx = range(3) # 0 -> x, 1 -> y, 2 -> z
        pp_st = H.makeSoltab(outSolsetPP, 'piercepoint', outSoltabPP,
            axesNames=['time', 'piercepoint', 'coord'], axesVals=[times_out, pp_indx,
            pp_pos_indx], vals=pp, weights=np.ones_like(pp))

        # Add histories
        sw = solWriter(tec_fit_st)
        sw.addHistory('CREATE (by TECSCREEN operation)')
        sw = solWriter(pp_st)
        sw.addHistory('CREATE (by TECSCREEN operation)')
        indx += 1

    return 0


