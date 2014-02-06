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

    N_sources = source_positions.shape[0]
    N_stations = station_positions.shape[0]
    N_piercepoints = N_stations * N_sources
    N_times = len(times)

    me = pyrap.measures.measures()
    position = me.position('ITRF', '%fm' % station_positions[0,0], '%fm' % station_positions[0,1], '%fm' % station_positions[0,2])
    me.doframe(position)


    pp = numpy.zeros((N_times, N_piercepoints,3))
    airmass = numpy.zeros((N_times, N_piercepoints))

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
    Fits a screen to given TEC values and save results to a new solset
    """

    N_stations = len(station_names)
    N_sources = len(source_names)
    N_times = len(times)

    tec_fit_all = zeros((N_times, N_sources, N_stations))

    A = concatenate([kron(eye(N_sources), ones((N_stations,1))),
        kron(ones((N_sources,1)), eye(N_stations))], axis=1)

    N_piercepoints = N_sources*N_stations
    P = eye(N_piercepoints) - dot(dot(A, pinv(dot(A.T, A))), A.T)
    N_freqs = 1

   # Initialize arrays
    parms = {}
    v = {}
    for station_name in station_names:
        for source_name in source_names:

            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = zeros((N_times, N_freqs), dtype=double)
            parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

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

        tec_fit = dot(U[:,:order], dot(pinvB, rr1))
        tec_fit_white = dot(pinv(C), tec_fit)
        tec_fit_all[k,:,:] = tec_fit_white.reshape((N_sources, N_stations))

    return tec_fit_all


def run( step, parset, H ):
    """
    Fits a screen to TEC values derived by the TECFIT operation.

    The TEC values must be stored in the 'tecfit' soltab of the specified
    solset(s).

    The results of the fit are stored in the 'piercepoints' and 'tecfitwhite'
    soltabs of the specified solset(s).
    """
    import numpy as np
    from h5parm import solFetcher, solWriter

    solsets = getParSolsets( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    # Load TEC values from TECFIT operations
    for solset in solsets:
        soltab = H.getSoltab(solset, 'tecfit')
        t = solFetcher(soltab)
        station_dict = H.getAnt(solset)
        station_names = station_dict.keys()
        station_positions = station_dict.values()
        source_dict = H.getAnt(solset)
        source_positions = source_dict.values()
        r, axis_vals = t.getValues()
        freqs = axis_vals['freqs']
        times = axis_vals['times']

        N_sources = len(source_names)
        N_times = len(times)
        N_stations = len(station_names)
        N_piercepoints = N_sources * N_stations
        rr = reshape(r.transpose([0,2,1]), [ N_piercepoints, N_times])

        # Find pierce points and airmass values for given screen height
        pp, airmass = calculate_piercepoints(station_positions, source_positions, times, height)

        # Fit a TEC screen
        r_0 = 10e3
        beta = 5.0/3.0
        tec_fit = fit_screen_to_tec(station_names, source_names, pp, airmass, rr, times, height, order, r_0, beta)

        # Write the results to the output solset
        dirs_out = {}
        for s in source_selection:
            dirs_out[source_names[s]] = source_positions[s]
        times_out = times
        ants_out = {}
        for s in station_selection:
            ants_out[station_names[s]] = station_positions[s]
        H.makeSoltab(solset, 'piercepoints' 'piercepoints', axesNames=['dir', 'time', 'ant'], \
                        axesVals=[dirs_out, times_out, ants_out], vals=pp)
        H.makeSoltab(solset, 'tecfitwhite' 'tecfitwhite', axesNames=['dir', 'time', 'ant'], \
                        axesVals=[dirs_out, times_out, ants_out], vals=tec_fit)

        # Store beta and r_0 as attributes of the tecfitwhite soltab
        tec_fit_st = H.getSoltab(solset, 'tecfitwhite')
        tec_fit_st.attrs['beta'] = beta
        tec_fit_st.attrs['r0'] = r_0

        # Add history
        sw = solWriter(tec_fit_st)
        sw.addHistory('CREATE (by TECSCREEN operation)')
        pp_st = H.getSoltab(solset, 'piercepoints')
        sw = solWriter(pp_st)
        sw.addHistory('CREATE (by TECSCREEN operation)')

    return 0


