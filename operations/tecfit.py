#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the phase-gradient TEC-fitting operation for LoSoTo

import logging
from operations_lib import *
from pylab import *
import lofar.expion
import pyrap.measures
import scipy.interpolate
import re

logging.debug('Loading TECFIT module.')


def collect_solutions(H, freq_tol=1e6):
    """
    Collects and returns phase solutions, etc. needed for fitting
    """
    logging.info("--> Collecting solutions needed for TEC fitting...")

    # Determine axis lengths
    sources = []
    freqs = []
    stations = []
    solsets = H.getSolsets().keys()
    N_times_max = 0
    for solset in solsets[:]:
        has_phase_st = False
        soltabs = H.getSoltabs(solset)

        for soltab in soltabs:
            # Use only soltabs of type phase
            if soltabs[soltab]._v_title == 'phase':
                has_phase_st = True
                solset_sources = H.getSou(solset)
                # Ignore direction-independent soltabs
                if len(solset_sources) == 1 and solset_sources[0] == 'pointing':
                    dir_indep = True
                else:
                    dir_indep = False
                if not dir_indep:
                    sources += H.getSou(solset)
                    stations += H.getAnt(solset)
                    t = solFetcher(soltabs[soltab])
                    solns, axis_vals = t.getValues()
                    times = axis_vals['time']
                    N_times = len(times)
                    freqs.append(axis_vals['freq'][0])
                    if N_times > N_times_max:
                        N_times_max = N_times
                        times_max = times
        if not has_phase_st:
            solsets.remove(solset)

    # Consolidate frequencies into bands (defined by freq_tol parameter)
    has_dup = True
    while has_dup:
        has_dup = False
        for freq in freqs[:]:
            nfreq = len( find( (array(freqs) > (freq-freq_tol)) &
                               (array(freqs) < (freq+freq_tol)) ) )
            if nfreq > 1:
                has_dup = True
                freqs.remove(freq)
                break
    freqs = array(sorted(freqs))
    N_freqs = len(freqs)
    sources_set = set(sources)
    sources_set.remove('pointing')
    source_names = array(sorted(list(sources_set)))
    N_sources = len(set(sources))
    N_stations = len(set(stations))
    N_times = N_times_max
    logging.info('Number of sources: {0}'.format(N_sources))
    logging.info('Number of stations: {0}'.format(N_stations))
    logging.info('Number of times: {0}'.format(N_times))
    logging.info('Number of freqs: {0}'.format(N_freqs))

    # Initialize the arrays
    freqwidths = zeros(N_freqs)
    timewidths = zeros(N_times)
    m = zeros((N_sources, N_freqs))
    phases0 = zeros((N_sources, N_stations, N_freqs, N_times))
    phases1 = zeros((N_sources, N_stations, N_freqs, N_times))
    flags = ones((N_sources, N_stations, N_freqs, N_times))
    source_positions = zeros((N_sources, 2))
    solset_array_dir_dep = chararray((N_sources, N_freqs), itemsize=100)
    solset_array_dir_indep = chararray((N_freqs), itemsize=100)

    # Populate the arrays
    for solset in solsets:
        source_dict = H.getSou(solset)
        sources = source_dict.keys()
        soltabs = H.getSoltabs(solset)

        for soltab in soltabs:
            t = solFetcher(soltabs[soltab])
            solns, axes = t.getValues()
            freq = axes['freq'][0]
            freq_indx = find( (array(freqs) > (freq+freq_tol)) & (array(freqs) < (freq-freq_tol)) )

            for source in sources:
                if source in source_names:
                    source_indx = find(source_names == source)
                    m[source_indx, freq_indx] = 1
                    solset_array_dir_dep[source_indx, freq_indx] = solset
                else:
                    solset_array_dir_indep[freq_indx] = solset

    # Collect station properties and pointing direction
    station_dict = H.getAnt(solsets[0])
    station_names = station_dict.keys()
    station_positions = np.zeros((len(station_names), 3), dtype=np.float)
    for i, station_name in enumerate(station_names):
        station_positions[i, 0] = station_dict[station_name][0]
        station_positions[i, 1] = station_dict[station_name][1]
        station_positions[i, 2] = station_dict[station_name][2]

    source_dict = H.getSou(solsets[0])
    ra = source_dict['pointing'][0]
    dec = source_dict['pointing'][1]
    pointing = array([cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)])

    # Collect source positions and phase solutions for each band
    for i, source1 in enumerate(source_names):
        for k in range(N_freqs):
            if m[i, k]:
                solset_name = solset_array_dir_dep[i, k]
                soltab = H.getSoltab(solset=solset_name, soltab='phase000')
                solset_dir_indep_name = solset_array_dir_indep[k]
                source_dict = H.getSou(solset_name)
                source_positions[i, ] = source_dict[source1]

                if solset_dir_indep_name != '':
                    soltab_dir_indep = H.getSoltab(solset=solset_dir_indep_name, soltab='phase000')
                    sf_dir_indep = solFetcher(soltab_dir_indep)
                else:
                    sf_dir_indep = None
                sf_dir_dep = solFetcher(soltab)

                for l, station in enumerate(station_names):
                    sf_dir_dep.setSelection(ant=station, pol='XX', dir=source1)
                    v1_phase, times_dir_dep = sf_dir_dep.getValues()
                    if sf_dir_indep is not None:
                        sf_dir_indep.setSelection(ant=station, pol='XX', dir=source1)
                        v1_dir_indep, times_dir_indep = sf_dir_indep.getValues()
                        v1_dir_indep_interp = interpolate_phase(v1_dir_indep, times_dir_indep, times_dir_dep)
                        v1_phase += v1_dir_indep_interp
                    if len(times_dir_dep) != N_times_max:
                        phases0 = interpolate_phase(phases0, times_dir_indep, times_max)
                    phases0[i, l, k, :] = v1_phase

                    sf_dir_dep.setSelection(ant=station, pol='YY', dir=source1)
                    v1_phase, times_dir_dep = sf_dir_dep.getValues()
                    if sf_dir_indep is not None:
                        sf_dir_indep.setSelection(ant=station, pol='YY', dir=source1)
                        v1_dir_indep, times_dir_indep = sf_dir_indep.getValues()
                        v1_dir_indep_interp = interpolate_phase(v1_dir_indep, times_dir_indep, times_dir_dep)
                        v1_phase += v1_dir_indep_interp
                    if len(times_dir_dep) != N_times_max:
                        phases1 = interpolate_phase(phases1, times_dir_indep, times_max)
                    phases1[i, l, k, :] = v1_phase

                    flags[i, l, k, :] = sf_dir_dep.getValues(weights=True, retAxesVals=False)

    return phases0, phases1, flags, m, station_names, station_positions, source_names, source_positions, freqs, times, pointing


def interpolate_phase(phase1, time1, time2, interpMethod='cubic'):
    """Interpolates phase solutions from one time grid to another

    If time2 is a coarser grid, averaging is done.
    """

    phase1 = unwrap_fft(phase1)

    if len(time2) < len(time1):
        # Average
        tratio = float(len(time2)) / float(len(time1))
        nslots = np.ceil(tratio)
        phase1, time1 = average_phase(phase1, nslots)

    # Interpolate
    valsNew = scipy.interpolate.griddata(time1, phase1, time2, interpMethod)

    return valsNew


def average_phase(phase, nslots, times=None):
    """Averages input phases over nslots time slots"""
    phase = (phase+numpy.pi) % (2*numpy.pi) - numpy.pi
    phase_avg = []
    p_subarray = []
    time_avg = []
    t_subarray = []
    for i, p in enumerate(phase):
        if i % nslots == 0 and len(p_subarray) > 0:
            phase_avg.append(numpy.mean(numpy.array(p_subarray)))
            p_subarray = []
        p_subarray.append(p)
        if i == len(phase)-1:
            phase_avg.append(numpy.mean(numpy.array(p_subarray)))

    if times is not None:
        time_avg = []
        t_subarray = []
        for i, t in enumerate(times):
            if i % nslots == 0 and len(t_subarray) > 0:
                time_avg.append(numpy.mean(numpy.array(t_subarray)))
                t_subarray = []
            t_subarray.append(t)
            if i == len(time)-1:
                time_avg.append(numpy.mean(numpy.array(t_subarray)))
        return numpy.array(phase_avg), numpy.array(time_avg)
    else:
        return numpy.array(phase_avg)


def fit_tec_per_source_pair(phases, flags, mask, freqs, init_sols = None, init_sols_per_pair = False, propagate = False, nband_min=2):

    sols_list = []
    eq_list = []

    N_sources = phases.shape[0]
    N_stations = phases.shape[1]
    N_times = phases.shape[3]

    if init_sols is None:
        init_sols = zeros((N_sources, N_times, N_stations), dtype = float)

    source_pairs = []

    k = 0
    for i in range(N_sources):
        for j in range(i):
            subband_selection = find(mask[i,:] * mask[j,:])
            if len(subband_selection) < nband_min:
                continue


            print i,j
            source_pairs.append((i,j))
            p = phases[i, :, subband_selection, :] - phases[j, :, subband_selection, :]
            A = zeros((len(subband_selection), 1))
            A[:,0] = 8.44797245e9/freqs[subband_selection]

            flags_source_pair = flags[i, :, subband_selection, :] * flags[j, :, subband_selection, :]
            constant_parms = zeros((1, N_stations), dtype = bool)
            sols = zeros((N_times, N_stations), dtype = float)

            if init_sols_per_pair:
                p_0 = init_sols[k, 0, :][newaxis,:]
            else:
                p_0 = (init_sols[i, 0, :] - init_sols[j, 0, :])[newaxis,:]

            for t_idx in range(N_times):
                x = p[:,:,t_idx].copy()
                f = flags_source_pair[:,:,t_idx].copy()
                if not propagate:
                    if init_sols_per_pair:
                        p_0 = init_sols[k, t_idx, :][newaxis,:]
                    else:
                        p_0 = (init_sols[i, t_idx, :] - init_sols[j, 0, :])[newaxis,:]
                sol = lofar.expion.baselinefitting.fit(x, A, p_0, f, constant_parms)
                if propagate:
                    p_0 = sol.copy()
                sols[t_idx, :] = sol[0,:]
            #subplot(1,2,2)
            sols = sols[:, :]-mean(sols[:, :], axis=1)[:,newaxis]
            weight = len(subband_selection)
            sols_list.append(sols*weight)
            eq = zeros(N_sources)
            eq[i] = weight
            eq[j] = -weight
            eq_list.append(eq)
            k += 1

    sols = array(sols_list)
    B = array(eq_list)
    source_selection = find(sum(abs(B), axis=0))
    N_sources = len(source_selection)
    if N_sources == 0:
        logging.error('No sources with enough bands.')
        return None, None, None, None
    pinvB = pinv(B)
    r = dot(pinvB, sols.transpose([1,0,2]))
    r = r[source_selection,:,:]

    return r, source_selection, sols, source_pairs


def add_stations(station_selection, phases0, phases1, flags, mask, station_names, station_positions, source_names, source_selection, times, freqs, r, nband_min=2):

    #No good solutions for RS409 RS310 RS208
    bad_stations = ['RS409LBA', 'RS310LBA', 'RS208LBA']

    N_sources_selected = len(source_selection)
    N_stations_selected = len(station_selection)
    N_piercepoints = N_sources_selected * N_stations_selected
    N_times = len(times)
    N_stations = len(station_names)
    N_sources = len(source_names)

    D = resize( station_positions, ( N_stations, N_stations, 3 ) )
    D = transpose( D, ( 1, 0, 2 ) ) - D
    D = sqrt(sum( D**2, axis=2 ))

    station_selection1 = station_selection

    stations_to_add = array([i for i in range(len(station_names)) if i not in station_selection1 and station_names[i] not in bad_stations])
    print station_names[stations_to_add]

    q = r

    while len(stations_to_add)>0:
        D1 = D[stations_to_add[:,newaxis], station_selection1[newaxis,:]]

        minimum_distance = amin(D1, axis=1)
        station_to_add = stations_to_add[argmin(minimum_distance)]
        print station_names[station_to_add]
        station_selection1 = append(station_selection1, station_to_add)
        N_stations_selected1 = len(station_selection1)

        # Remove station from list
        stations_to_add = stations_to_add[stations_to_add != station_to_add]

        sols_list = []
        eq_list = []
        min_e_list = []

        for ii, i in enumerate(source_selection):
            for jj, j in enumerate(source_selection):
                if j == i:
                    break
                subband_selection = find(mask[i,:] * mask[j,:])
                if len(subband_selection) < nband_min:
                    continue
                print '===================='
                print i,j, 'number of subbands = ', len(subband_selection)
                print '===================='
                p0 = phases0[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] - phases0[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                p0 = p0 - mean(p0, axis=0)[newaxis,:,:]
                p1 = phases1[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] - phases1[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                p1 = p1 - mean(p1, axis=0)[newaxis,:,:]
                A = zeros((len(subband_selection), 1))
                A[:,0] = 8.44797245e9/freqs[subband_selection]

                flags_source_pair = flags[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] * flags[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                constant_parms = zeros((1, N_stations_selected1), dtype = bool)
                sols = zeros((N_times, N_stations_selected1), dtype = float)

                for t_idx in range(N_times):
                    min_e = Inf
                    for offset in linspace(-0.1, 0.1,21) :
                        p_0 = zeros((1, N_stations_selected1), double)
                        p_0[0,:N_stations_selected1-1] = (q[ii, t_idx, :] - q[jj, t_idx, :])[newaxis,:]
                        p_0[0,-1] = offset

                        x = p0[:,:,t_idx].copy()
                        f = flags_source_pair[:,:,t_idx].copy()
                        sol0 = lofar.expion.baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                        sol0 -= mean(sol0)
                        residual = mod(dot(A,sol0) - x.T + pi, 2*pi) - pi
                        residual = residual[f.T==0]
                        e = var(residual)

                        x = p1[:,:,t_idx].copy()
                        f = flags_source_pair[:,:,t_idx].copy()
                        sol1 = lofar.expion.baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                        sol1 -= mean(sol1)
                        residual = mod(dot(A,sol1) - x.T + pi, 2*pi) - pi
                        residual = residual[f.T==0]
                        e += var(residual)


                        if e<min_e:
                            min_e = e
                            sols[t_idx, :] = (sol0[0,:] + sol1[0,:])/2

                        #residual = mod(dot(A,sol) - x.T + pi, 2*pi) - pi
                        #residual = residual[f.T==0]
                        #ee = var(residual)
                        #e += ee
                        #e_list.append(sum(abs(residual) > 3*std(residual)))
                    #print offset, e

                #weight = len(subband_selection) / min_e

                ### Remove outliers
                for kk in range(10):
                    s = sols[:,-1].copy()
                    selection = zeros(len(s), bool)
                    for t_idx in range(len(s)):
                        start_idx = max(t_idx-10, 0)
                        end_idx = min(t_idx+10, len(s))
                        selection[t_idx] = sum(abs(s[start_idx:end_idx] - s[t_idx])<0.02) > (end_idx - start_idx - 8)
                    outliers = find(logical_not(selection))
                    if len(outliers) == 0:
                        break
                    for t_idx in outliers:
                        print "outlier", t_idx
                        try:
                            idx0 = find(selection[:t_idx])[-1]
                        except IndexError:
                            idx0 = -1
                        try:
                            idx1 = find(selection[t_idx+1:])[0]+t_idx+1
                        except IndexError:
                            idx1 = -1
                        if idx0 == -1:
                            s[t_idx] = s[idx1]
                        elif idx1 == -1:
                            s[t_idx] = s[idx0]
                        else:
                            s[t_idx] = (s[idx0] * (idx1-t_idx) + s[idx1] * (t_idx-idx0)) / (idx1-idx0)

                        p_0 = zeros((1, N_stations_selected1), double)
                        p_0[0,:] = sols[t_idx,:]
                        p_0[0,-1] = s[t_idx]

                        x = p0[:,:,t_idx].copy()
                        f = flags_source_pair[:,:,t_idx].copy()
                        sol0 = lofar.expion.baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                        sol0 -= mean(sol0)

                        x = p1[:,:,t_idx].copy()
                        sol1 = lofar.expion.baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                        sol1 -= mean(sol1)

                        print sols[t_idx, -1], (sol0[0,-1] + sol1[0,-1])/2, s[t_idx]
                        sols[t_idx, :] = (sol0[0,:] + sol1[0,:])/2
                        #sols[t_idx, -1] = s[t_idx]

                weight = 1.0
                sols_list.append(weight*sols)
                min_e_list.append(min_e)
                eq = zeros(N_sources)
                eq[ii] = weight
                eq[jj] = -weight
                eq_list.append(eq)

        sols = array(sols_list)
        B = array(eq_list)
        pinvB = pinv(B)

        q = dot(pinvB, sols.transpose([1,0,2]))
        if N_stations_selected1 == 28:
            break

    return station_selection1, sols, q


def run( step, parset, H ):
    """
    Fit phase gradients to obtain TEC values per direction.

    Phase solutions are assumed to be stored in solsets of the H5parm file, one
    solset per band per field.

    The TEC values are stored in the 'tecfit' soltab of the 'ionosphere'
    solset.

    """
    import numpy as np
    from h5parm import solFetcher, solWriter

    solsets = getParSolsets( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )
    nband_min = getParAxis( step, parset, H, 'MinBands' )
    outSolset = 'ionosphere'
    H.makeSolset(outSolset)

    # Collect solutions, etc. into arrays for fitting.
    # At the moment, all solsets of H are assumed to be needed
    phases0, phases1, flags, mask, station_names, station_positions, source_names, source_positions, freqs, times, pointing = collect_solutions(H)

    # Build list of stations to include
    included_stations = []
    for ant in ants:
        included_stations += [s for s in station_names if re.search(ant, s)]
    excluded_stations = [s for s in station_names if s not in included_stations]

    # Select stations to use for first pass
    logging.info("Excluding stations: {0}".format(sort(excluded_stations)))
    mean_position = np.array([median(station_positions[:, 0]), median(station_positions[:, 1]), median(station_positions[:, 2])])
    station_selection1 = find(sqrt(sum((station_positions - mean_position)**2, axis=1)) < 2e3)
    station_selection = array([i for i in range(len(station_names)) if i in station_selection1 and station_names[i] not in excluded_stations])

    # Fit a TEC value to the phase solutions per source pair
    r0, source_selection, sols0, source_pairs = fit_tec_per_source_pair(phases0[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True, nband_min=nband_min)
    r1, source_selection, sols1, source_pairs = fit_tec_per_source_pair(phases1[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True, nband_min=nband_min)

    if r0 is None or r1 is None:
        return 1

    # take the mean of the two polarizations
    r = (r0+r1)/2

    # Try to add more stations
    station_selection, sols, r = add_stations(station_selection, phases0, phases1, flags, mask, station_names, station_positions, source_names, source_selection, times, freqs, r, nband_min=nband_min)

    # Save TEC values to the output solset
    dirs_out = {}
    for s in source_selection:
        dirs_out[source_names[s]] = source_positions[s]
    times_out = times
    ants_out = {}
    for s in station_selection:
        ants_out[station_names[s]] = station_positions[s]
    H.makeSoltab(outSolset, 'tecfit', 'tec', axesNames=['dir', 'time', 'ant'], \
                    axesVals=[dirs_out, times_out, ants_out], vals=r)

    # Add history
    tf_st = H.getSoltab(solset, 'tecfit')
    sw = solWriter(tf_st)
    sw.addHistory('CREATE (by TECFIT operation)')

    return 0


