#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the phase-gradient TEC-fitting operation for LoSoTo

import logging
from operations_lib import *
from pylab import *
from lofar.expion import baselinefitting
import pyrap.measures
import scipy.interpolate
import re

logging.debug('Loading TECFIT module.')


def collect_solutions(H, dirs=None, freq_tol=1e6):
    """
    Collects and returns phase solutions, etc. needed for fitting
    """
    logging.info("Scanning for solutions needed for TEC fitting...")

    # Determine axis lengths
    sources = []
    freqs = []
    stations = []
    solsets = H.getSolsets().keys()
    N_times_max = 0
    first_solset = None
    for solset in solsets[:]:
        if 'ion' in solset:
            continue
        logging.info('  --Scanning solution set {0}...'.format(solset))
        has_phase_st = False
        soltabs = H.getSoltabs(solset)

        for soltab in soltabs:
            # Only use soltabs with phase solutions
            if 'phase' in soltabs[soltab]._v_title:
                logging.info('    --Scanning solution table {0}...'.format(soltab))
                has_phase_st = True
                solset_sources = H.getSou(solset)
                # Ignore direction-independent soltabs
                if len(solset_sources) == 1 and 'pointing' in solset_sources:
                    logging.info('      Found direction-independent solutions')
                    dir_indep = True
                    soln_type_dirindep = soltabs[soltab]._v_title
                    if first_solset is None:
                        first_solset = solset
                else:
                    logging.info('      Found direction-dependent solutions')
                    dir_indep = False
                if not dir_indep:
                    soln_type_dirdep = soltabs[soltab]._v_title
                    logging.info('      Found sources: {0}'.format(H.getSou(solset).keys()))
                    sources += H.getSou(solset).keys()
                    stations += H.getAnt(solset).keys()
                    t = solFetcher(soltabs[soltab])
                    solns, axis_vals = t.getValues()
                    times = axis_vals['time']
                    N_times = len(times)
                    freqs.append(axis_vals['freq'][0])
                    if N_times > N_times_max:
                        N_times_max = N_times
                        times_max = times
        if not has_phase_st:
            logging.info('    --No phase solutions found')
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
    stations_set = set(stations)
    sources_set = set(sources)
    if 'pointing' in sources_set:
        sources_set.remove('pointing')
    if dirs is not None:
        sources_filt = []
        for dir in dirs:
            if dir in sources_set:
                sources_filt.append(dir)
        sources_set = set(sources_filt)
    source_names = array(sorted(list(sources_set)))
    N_sources = len(set(source_names))
    N_stations = len(set(stations))
    N_times = N_times_max
    logging.info('Scanning complete')
    logging.info('  Number of sources: {0}'.format(N_sources))
    logging.info('  Number of stations: {0}'.format(N_stations))
    logging.info('  Number of times: {0}'.format(N_times))
    logging.info('  Number of freqs: {0}'.format(N_freqs))

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
        if 'ion' in solset:
            continue
        source_dict = H.getSou(solset)
        sources = source_dict.keys()
        soltabs = H.getSoltabs(solset)
        # TODO: restrict to soltabs with phase and amp only (e.g., exclude tec soltabs)

        for soltab in soltabs:
            t = solFetcher(soltabs[soltab])
            solns, axes = t.getValues()
            freq = axes['freq'][0]
            freq_indx = find( (array(freqs) > (freq-freq_tol)) & (array(freqs) < (freq+freq_tol)) )

            for source in sources:
                if source in source_names:
                    source_indx = find(source_names == source)
                    m[source_indx, freq_indx] = 1
                    solset_array_dir_dep[source_indx, freq_indx] = solset
                elif source=='pointing' and len(sources) == 1:
                    solset_array_dir_indep[freq_indx] = solset

    # Collect station properties and pointing direction
    station_dict = H.getAnt(first_solset)
    station_names = np.array(list(stations_set))
    station_positions = np.zeros((len(station_names), 3), dtype=np.float)
    for i, station_name in enumerate(station_names):
        station_positions[i, 0] = station_dict[station_name][0]
        station_positions[i, 1] = station_dict[station_name][1]
        station_positions[i, 2] = station_dict[station_name][2]

    source_dict = H.getSou(first_solset)
    ra = source_dict['pointing'][0]
    dec = source_dict['pointing'][1]
    pointing = array([cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)])

    # Collect source positions and phase solutions for each band
    logging.info('Collecting phase solutions...')
    for i, source1 in enumerate(source_names):
        logging.info('  Source {0}'.format(source1))
        for k in range(N_freqs):
            if m[i, k]:
                solset_name = solset_array_dir_dep[i, k]
                soltab = H.getSoltab(solset=solset_name, soltab=soln_type_dirdep+'000')
                solset_dir_indep_name = solset_array_dir_indep[k]
                source_dict = H.getSou(solset_name)
                source_positions[i, ] = source_dict[source1]

                if solset_dir_indep_name != '':
                    soltab_dir_indep = H.getSoltab(solset=solset_dir_indep_name, soltab=soln_type_dirindep+'000')
                    sf_dir_indep = solFetcher(soltab_dir_indep)
                else:
                    sf_dir_indep = None
                sf_dir_dep = solFetcher(soltab)

                for l, station in enumerate(station_names):
                    if soln_type_dirdep == 'scalarphase':
                        sf_dir_dep.setSelection(ant=station, dir=source1)
                    else:
                        sf_dir_dep.setSelection(ant=station, pol='XX', dir=source1)
                    values_dir_dep = sf_dir_dep.getValues()
                    v1_phase = array(values_dir_dep[0]).squeeze()
                    times_dir_dep = values_dir_dep[1]['time']
                    if sf_dir_indep is not None:
                        if soln_type_dirindep == 'scalarphase':
                            sf_dir_indep.setSelection(ant=station)
                        else:
                            sf_dir_indep.setSelection(ant=station, pol='XX')
                        values_dir_indep = sf_dir_indep.getValues()
                        v1_dir_indep = array(values_dir_indep[0]).squeeze()
                        times_dir_indep = values_dir_indep[1]['time']
                        v1_dir_indep_interp = interpolate_phase(v1_dir_indep, times_dir_indep, times_dir_dep)
                        v1_phase += v1_dir_indep_interp
                    if len(times_dir_dep) != N_times_max:
                        phases0 = interpolate_phase(phases0, times_dir_indep, times_max)
                    phases0[i, l, k, :] = v1_phase

                    if soln_type_dirdep == 'scalarphase':
                        phases1[i, l, k, :] = v1_phase
                    else:
                        sf_dir_dep.setSelection(ant=station, pol='YY', dir=source1)
                        values_dir_dep = sf_dir_dep.getValues()
                        v1_phase = array(values_dir_dep[0]).squeeze()
                        times_dir_dep = values_dir_dep[1]['time']
                        if sf_dir_indep is not None:
                            if soln_type_dirindep == 'scalarphase':
                                sf_dir_indep.setSelection(ant=station)
                            else:
                                sf_dir_indep.setSelection(ant=station, pol='YY')
                            values_dir_indep = sf_dir_indep.getValues()
                            v1_dir_indep = array(values_dir_indep[0]).squeeze()
                            times_dir_indep = values_dir_indep[1]['time']
                            v1_dir_indep_interp = interpolate_phase(v1_dir_indep, times_dir_indep, times_dir_dep)
                            v1_phase += v1_dir_indep_interp
                        if len(times_dir_dep) != N_times_max:
                            phases1 = interpolate_phase(phases1, times_dir_indep, times_max)
                        phases1[i, l, k, :] = v1_phase

                    flags[i, l, k, :] = sf_dir_dep.getValues(weight=True, retAxesVals=False)

    # Invert the weights to give flags (0 => unflagged, 1 => flagged)
    zeroflags = np.where(flags == 0.0)
    oneflags = flags.nonzero()
    flags[zeroflags] = 1.0
    flags[oneflags] = 0.0
    return (phases0, phases1, flags, m, station_names, station_positions,
        source_names, source_positions, freqs, times, pointing, soln_type_dirdep)


def interpolate_phase(phase1, time1, time2, interpMethod='cubic'):
    """Interpolates phase solutions from one time grid to another

    If time2 is a coarser grid, averaging is done instead.
    """

    phase1 = unwrap_fft(phase1)

    if len(time2) < len(time1):
        # Average
        tratio = float(len(time1)) / float(len(time2))
        nslots = np.ceil(tratio)
        valsNew = average_phase(phase1, nslots)
    else:
        # Interpolate
        valsNew = scipy.interpolate.griddata(time1, phase1, time2, interpMethod)

    return valsNew


def unwrap_fft(phase, iterations=1):
    """
    Unwrap phase using Fourier techniques.

    For details, see:
    Marvin A. Schofield & Yimei Zhu, Optics Letters, 28, 14 (2003)
    """

    puRadius=lambda x : np.roll( np.roll(
          np.add.outer( np.arange(-x.shape[0]/2+1,x.shape[0]/2+1)**2.0,
                        np.arange(-x.shape[1]/2+1,x.shape[1]/2+1)**2.0 ),
          x.shape[1]/2+1,axis=1), x.shape[0]/2+1,axis=0)+1e-9

    idt,dt=np.fft.ifft2,np.fft.fft2
    puOp=lambda x : idt( np.where(puRadius(x)==1e-9,1,puRadius(x)**-1.0)*dt(
          np.cos(x)*idt(puRadius(x)*dt(np.sin(x)))
         -np.sin(x)*idt(puRadius(x)*dt(np.cos(x))) ) )

    def phaseUnwrapper(ip):
       mirrored=np.zeros([x*2 for x in ip.shape])
       mirrored[:ip.shape[0],:ip.shape[1]]=ip
       mirrored[ip.shape[0]:,:ip.shape[1]]=ip[::-1,:]
       mirrored[ip.shape[0]:,ip.shape[1]:]=ip[::-1,::-1]
       mirrored[:ip.shape[0],ip.shape[1]:]=ip[:,::-1]

       return (ip+2*np.pi*
             np.round((puOp(mirrored).real[:ip.shape[0],:ip.shape[1]]-ip)
             /2/np.pi))

    phase2D = phase[:, None]
    i = 0
    if iterations < 1:
        interations = 1
    while i < iterations:
        i += 1
        phase2D = phaseUnwrapper(phase2D)

    return phase2D[:, 0]


def average_phase(phase, nslots, times=None):
    """Averages input phases over nslots time slots"""
    phase = (phase+np.pi) % (2*np.pi) - np.pi
    phase_avg = []
    p_subarray = []
    time_avg = []
    t_subarray = []
    for i, p in enumerate(phase):
        if i % nslots == 0 and len(p_subarray) > 0:
            phase_avg.append(np.mean(np.array(p_subarray)))
            p_subarray = []
        p_subarray.append(p)
        if i == len(phase)-1:
            phase_avg.append(np.mean(np.array(p_subarray)))

    if times is not None:
        time_avg = []
        t_subarray = []
        for i, t in enumerate(times):
            if i % nslots == 0 and len(t_subarray) > 0:
                time_avg.append(np.mean(np.array(t_subarray)))
                t_subarray = []
            t_subarray.append(t)
            if i == len(time)-1:
                time_avg.append(np.mean(np.array(t_subarray)))
        return np.array(phase_avg), np.array(time_avg)
    else:
        return np.array(phase_avg)


def fit_tec_per_source_pair(phases, flags, mask, freqs, init_sols=None,
    init_sols_per_pair=False, propagate=False, nband_min=2, offset=0):

    sols_list = []
    eq_list = []
    vars_list = []
    var_eq_list = []

    N_sources = phases.shape[0]
    N_stations = phases.shape[1]
    N_times = phases.shape[3]
    k = 0
    for i in range(N_sources):
        for j in range(i):
            k += 1
    N_pairs = k

    if init_sols is None and not init_sols_per_pair:
        init_sols = zeros((N_sources, N_times, N_stations), dtype = float) + offset
    elif init_sols is None and init_sols_per_pair:
        init_sols = zeros((N_pairs, N_times, N_stations), dtype = float) + offset
    vars = zeros((N_times, N_stations), dtype = float)
    source_pairs = []

    k = 0
    for i in range(N_sources):
        for j in range(i):
            subband_selection = find(mask[i,:] * mask[j,:])
            if len(subband_selection) < nband_min:
                continue
            logging.info('Fitting TEC values for source pair: {0}-{1}'.format(i, j))
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
                sol = baselinefitting.fit(x, A, p_0, f, constant_parms)

                if propagate:
                    p_0 = sol.copy()
                sols[t_idx, :] = sol[0, :]
            sols = sols[:, :]-mean(sols[:, :], axis=1)[:,newaxis]
            weight = len(subband_selection)
            sols_list.append(sols*weight)
            vars_list.append(vars*weight)
            eq = zeros(N_sources)
            eq[i] = weight
            eq[j] = -weight
            eq_list.append(eq)
            var_eq = zeros(N_sources)
            var_eq[i] = weight
            var_eq[j] = weight
            var_eq_list.append(var_eq)
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
    r = r[source_selection, :, :]

    return r, source_selection, sols, source_pairs


def add_stations(station_selection, phases0, phases1, flags, mask, station_names, station_positions, source_names, source_selection, times, freqs, r, nband_min=2, soln_type='phase', nstations_max=None):

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

    stations_to_add = array([i for i in range(len(station_names)) if i not in station_selection1])

    # Check if desired number of stations is already reached
    if nstations_max is not None:
        if len(station_selection1) >= nstations_max:
            return station_selection1, None, r

    q = r
    while len(stations_to_add)>0:
        D1 = D[stations_to_add[:,newaxis], station_selection1[newaxis,:]]

        minimum_distance = amin(D1, axis=1)
        station_to_add = stations_to_add[argmin(minimum_distance)]
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
                logging.info('Fitting {0} TEC values for source pair: {1}-{2}'.format(station_names[station_to_add], i, j))
                p0 = phases0[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] - phases0[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                p0 = p0 - mean(p0, axis=0)[newaxis,:,:]
                if soln_type != 'scalarphase':
                    p1 = phases1[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] - phases1[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                    p1 = p1 - mean(p1, axis=0)[newaxis,:,:]
                A = zeros((len(subband_selection), 1))
                A[:,0] = 8.44797245e9/freqs[subband_selection]

                flags_source_pair = flags[i, station_selection1[:,newaxis], subband_selection[newaxis,:], :] * flags[j, station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                constant_parms = zeros((1, N_stations_selected1), dtype = bool)
                sols = zeros((N_times, N_stations_selected1), dtype = float)

                for t_idx in range(N_times):
                    if np.mod(t_idx, 10) == 0:
                        min_e = Inf
                        for offset in linspace(-0.1, 0.1,21) :
                            p_0 = zeros((1, N_stations_selected1), double)
                            p_0[0,:N_stations_selected1-1] = (q[ii, t_idx, :] - q[jj, t_idx, :])[newaxis,:]
                            p_0[0,-1] = offset

                            x = p0[:,:,t_idx].copy()
                            f = flags_source_pair[:,:,t_idx].copy()
                            sol0 = baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                            sol0 -= mean(sol0)
                            residual = mod(dot(A,sol0) - x.T + pi, 2*pi) - pi
                            residual = residual[f.T==0]
                            e = var(residual)

                            if soln_type != 'scalarphase':
                                x = p1[:,:,t_idx].copy()
                                f = flags_source_pair[:,:,t_idx].copy()
                                sol1 = baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                                sol1 -= mean(sol1)
                                residual = mod(dot(A,sol1) - x.T + pi, 2*pi) - pi
                                residual = residual[f.T==0]
                                e += var(residual)
                            else:
                                sol1 = sol0

                            if e<min_e:
                                min_e = e
                                p_0_best = p_0
                                sols[t_idx, :] = (sol0[0,:] + sol1[0,:])/2
                    else:
                        # Use previous init
                        x = p0[:,:,t_idx].copy()
                        f = flags_source_pair[:,:,t_idx].copy()
                        sol0 = baselinefitting.fit(x.T, A, p_0_best, f, constant_parms)
                        sol0 -= mean(sol0)

                        if soln_type != 'scalarphase':
                            x = p1[:,:,t_idx].copy()
                            f = flags_source_pair[:,:,t_idx].copy()
                            sol1 = baselinefitting.fit(x.T, A, p_0_best, f, constant_parms)
                            sol1 -= mean(sol1)
                        else:
                            sol1 = sol0
                        sols[t_idx, :] = (sol0[0,:] + sol1[0,:])/2

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
                        sol0 = baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                        sol0 -= mean(sol0)

                        if soln_type != 'scalarphase':
                            x = p1[:,:,t_idx].copy()
                            sol1 = baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                            sol1 -= mean(sol1)
                        else:
                            sol1 = sol0

                        sols[t_idx, :] = (sol0[0,:] + sol1[0,:])/2

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
        if nstations_max is not None:
            if N_stations_selected1 == nstations_max:
                break

    return station_selection1, sols, q


def run( step, parset, H ):
    """
    Fit source-to-source phase gradients to obtain TEC values per direction.

    Phase solutions are assumed to be stored in solsets of the H5parm file, one
    solset per band per field.

    The TEC values are stored in the specified output soltab with type 'tec'.

    """
    import numpy as np
    from h5parm import solFetcher, solWriter

    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )
    nband_min = int(parset.getString('.'.join(["LoSoTo.Steps", step, "MinBands"]), '8' ))
    nstations_max = int(parset.getString('.'.join(["LoSoTo.Steps", step, "MaxStations"]), '100' ))
    outSolset = parset.getString('.'.join(["LoSoTo.Steps", step, "OutSoltab"]), '' ).split('/')[0]
    outSoltab = parset.getString('.'.join(["LoSoTo.Steps", step, "OutSoltab"]), '' ).split('/')[1]

    # Collect solutions, etc. into arrays for fitting.
    # At the moment, all solsets of H are searched for phase solutions.
    phases0, phases1, flags, mask, station_names, station_positions, source_names, source_positions, freqs, times, pointing, soln_type = collect_solutions(H, dirs=dirs)

    # Build list of stations to include
    included_stations = []
    if ants is not None:
        for ant in ants:
            included_stations += [s for s in station_names if re.search(ant, s)]
    else:
        included_stations = station_names
    excluded_stations = [s for s in station_names if s not in included_stations]

    # Select stations to use for first pass
    dist_cut_m = 2e3
    logging.info("Excluding stations: {0}".format(sort(excluded_stations)))
    mean_position = np.array([median(station_positions[:, 0]), median(station_positions[:, 1]), median(station_positions[:, 2])])
    station_selection1 = find(sqrt(sum((station_positions - mean_position)**2, axis=1)) < dist_cut_m)
    station_selection = array([i for i in range(len(station_names)) if i in station_selection1 and station_names[i] not in excluded_stations])

    # Fit a TEC value to the phase solutions per source pair. No search for the
    # global minimum is done
    search_tec = True
    nsteps = 21
    if soln_type == 'scalarphase':
        r, source_selection, sols0, source_pairs = fit_tec_per_source_pair(phases0[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True, nband_min=nband_min)
        if r is None:
            return 1
    else:
        r0, source_selection, sols0, source_pairs = fit_tec_per_source_pair(phases0[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True, nband_min=nband_min)
        r1, source_selection, sols1, source_pairs = fit_tec_per_source_pair(phases1[:,station_selection,:,:], flags[:,station_selection,:,:], mask, freqs, propagate = True, nband_min=nband_min)

        if r0 is None or r1 is None:
            return 1

        # take the mean of the two polarizations
        r = (r0+r1)/2

    # Add stations by searching for global minimum in solution space
    station_selection, sols, r = add_stations(station_selection, phases0, phases1, flags, mask, station_names, station_positions, source_names, source_selection, times, freqs, r, nband_min=nband_min, soln_type=soln_type, nstations_max=nstations_max)

    # Save TEC values to the output solset
    solset = H.makeSolset(outSolset)
    dirs_out = []
    dirs_pos = []
    for s in source_selection:
        dirs_out.append(source_names[s])
        dirs_pos.append(source_positions[s])
    sourceTable = solset._f_get_child('source')
    sourceTable.append(zip(*(dirs_out, dirs_pos)))

    times_out = times

    ants_out = []
    ants_pos = []
    for s in station_selection:
        ants_out.append(station_names[s])
        ants_pos.append(station_positions[s])
    antennaTable = solset._f_get_child('antenna')
    antennaTable.append(zip(*(ants_out, ants_pos)))

    tf_st = H.makeSoltab(outSolset, 'tec', outSoltab, axesNames=['dir', 'time', 'ant'],
        axesVals=[dirs_out, times, ants_out], vals=r, weights=np.ones_like(r))

    # Add history
    sw = solWriter(tf_st)
    sw.addHistory('CREATE (by TECFIT operation)')

    return 0


