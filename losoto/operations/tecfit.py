#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is the phase-gradient TEC-fitting operation for LoSoTo

import logging
from losoto.operations_lib import *

logging.debug('Loading TECFIT module.')


def collect_solutions(H, dirs=None, freq_tol=1e6, solsets=None):
    """
    Collects and returns phase solutions, etc. needed for fitting

    Keyword arguments:
    H -- H5parm object
    dirs -- list of directions to use (None => all)
    freq_tol -- tolerance for grouping of phase solutions into bands (Hz)
    solsets -- list of solution sets in H to search (None => all)
    """
    import numpy as np
    from pylab import find
    import re
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar
    logging.info("Scanning for solutions needed for TEC fitting...")

    # Determine axis lengths
    sources = []
    freqs = []
    stations = []
    if solsets is None:
        solsets = H.getSolsets().keys()
    N_times_max = 0
    first_solset = None
    for solset in solsets[:]:
        logging.info('  --Scanning solution set {0}...'.format(solset))
        has_phase_st = False
        soltabs = H.getSoltabs(solset)

        for soltab in soltabs:
            # Only use soltabs with phase solutions
            if 'phase' in soltabs[soltab]._v_title:
                logging.info('    --Scanning solution table {0}...'.format(soltab))
                has_phase_st = True
                solset_sources = H.getSou(solset)
                if len(solset_sources) == 1 and 'pointing' in solset_sources:
                    logging.info('      Found direction-independent solutions')
                    dir_indep = True
                    soln_type_dirindep = soltabs[soltab]._v_title
                    if first_solset is None:
                        first_solset = solset
                else:
                    logging.info('      Found direction-dependent solutions')
                    dir_indep = False
                    if first_solset is None:
                        first_solset = solset
                if not dir_indep:
                    soln_type_dirdep = soltabs[soltab]._v_title
                    logging.info('      Found sources: {0}'.format(
                        H.getSou(solset).keys()))
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
            nfreq = len( find( (np.array(freqs) > (freq-freq_tol)) &
                               (np.array(freqs) < (freq+freq_tol)) ) )
            if nfreq > 1:
                has_dup = True
                freqs.remove(freq)
                break
    freqs = np.array(sorted(freqs))
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
    source_names = np.array(sorted(list(sources_set)))
    N_sources = len(set(source_names))
    N_stations = len(set(stations))
    N_times = N_times_max
    times = times_max
    logging.info('Scanning complete')
    logging.info('  Number of sources: {0}'.format(N_sources))
    logging.info('  Number of stations: {0}'.format(N_stations))
    logging.info('  Number of times: {0}'.format(N_times))
    logging.info('  Number of bands: {0}'.format(N_freqs))
    if N_sources == 0 or N_stations == 0 or N_times == 0 or N_freqs == 0:
        logging.error('No solutions found.')
        return (None, None, None, None, None, None, None, None, None, None,
            None, None)

    # Initialize the arrays
    freqwidths = np.zeros(N_freqs)
    timewidths = np.zeros(N_times)
    m = np.zeros((N_sources, N_freqs))
    phases0 = np.zeros((N_sources, N_stations, N_freqs, N_times))
    phases1 = np.zeros((N_sources, N_stations, N_freqs, N_times))
    flags = np.ones((N_sources, N_stations, N_freqs, N_times))
    source_positions = np.zeros((N_sources, 2))
    solset_array_dir_dep = np.zeros((N_sources, N_freqs), dtype='|S100')
    solset_array_dir_indep = np.zeros(N_freqs, dtype='|S100')

    # Populate the arrays
    for solset in solsets:
        source_dict = H.getSou(solset)
        sources = source_dict.keys()
        soltabs = H.getSoltabs(solset)

        for soltab in soltabs:
            t = solFetcher(soltabs[soltab])
            solns, axes = t.getValues()
            freq = axes['freq'][0]
            freq_indx = find( (np.array(freqs) > (freq-freq_tol)) &
                (np.array(freqs) < (freq+freq_tol)) )

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
    pointing = np.array([np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec),
        np.sin(dec)])

    # Collect source positions and phase solutions for each band
    logging.info('Collecting phase solutions...')
    pbar = progressbar.ProgressBar(maxval=len(source_names)).start()
    for i, source1 in enumerate(source_names):
        for k in xrange(N_freqs):
            if m[i, k]:
                solset_name = str(solset_array_dir_dep[i, k])
                soltab = H.getSoltab(solset=solset_name, soltab=soln_type_dirdep+'000')
                solset_dir_indep_name = str(solset_array_dir_indep[k])
                source_dict = H.getSou(solset_name)
                source_positions[i, ] = source_dict[source1]

                if solset_dir_indep_name != '':
                    soltab_dir_indep = H.getSoltab(solset=solset_dir_indep_name,
                        soltab=soln_type_dirindep+'000')
                    sf_dir_indep = solFetcher(soltab_dir_indep)
                else:
                    sf_dir_indep = None
                sf_dir_dep = solFetcher(soltab)

                for l, station in enumerate(station_names):
                    if soln_type_dirdep == 'scalarphase':
                        sf_dir_dep.setSelection(ant=[station], dir=[source1])
                    else:
                        sf_dir_dep.setSelection(ant=[station], pol=['XX'], dir=[source1])
                    values_dir_dep = sf_dir_dep.getValues()
                    v1_phase = np.array(values_dir_dep[0]).squeeze()
                    times_dir_dep = values_dir_dep[1]['time']
                    ind = np.where(~np.isnan(v1_phase))
                    v1_phase = v1_phase[ind]
                    times_dir_dep = times_dir_dep[ind]

                    if sf_dir_indep is not None:
                        if soln_type_dirindep == 'scalarphase':
                            sf_dir_indep.setSelection(ant=[station])
                        else:
                            sf_dir_indep.setSelection(ant=[station], pol=['XX'])
                        values_dir_indep = sf_dir_indep.getValues()
                        v1_dir_indep = np.array(values_dir_indep[0]).squeeze()
                        times_dir_indep = values_dir_indep[1]['time']
                        ind = np.where(~np.isnan(v1_dir_indep))
                        v1_dir_indep = v1_dir_indep[ind]
                        times_dir_indep = times_dir_indep[ind]
                        v1_dir_indep_interp = interpolate_phase(v1_dir_indep,
                            times_dir_indep, times_dir_dep)
                        v1_phase += v1_dir_indep_interp
                    if len(times_dir_dep) != N_times_max:
                         v1_phase = interpolate_phase(v1_phase, times_dir_dep,
                            times_max)
                    phases0[i, l, k, :] = v1_phase

                    if soln_type_dirdep == 'scalarphase':
                        phases1[i, l, k, :] = v1_phase
                    else:
                        sf_dir_dep.setSelection(ant=[station], pol=['YY'], dir=[source1])
                        values_dir_dep = sf_dir_dep.getValues()
                        v1_phase = np.array(values_dir_dep[0]).squeeze()
                        times_dir_dep = values_dir_dep[1]['time']
                        ind = np.where(~np.isnan(v1_phase))
                        v1_phase = v1_phase[ind]
                        times_dir_dep = times_dir_dep[ind]

                        if sf_dir_indep is not None:
                            if soln_type_dirindep == 'scalarphase':
                                sf_dir_indep.setSelection(ant=[station])
                            else:
                                sf_dir_indep.setSelection(ant=[station], pol=['YY'])
                            values_dir_indep = sf_dir_indep.getValues()
                            v1_dir_indep = np.array(values_dir_indep[0]).squeeze()
                            times_dir_indep = values_dir_indep[1]['time']
                            ind = np.where(~np.isnan(v1_dir_indep))
                            v1_dir_indep = v1_dir_indep[ind]
                            times_dir_indep = times_dir_indep[ind]
                            v1_dir_indep_interp = interpolate_phase(v1_dir_indep,
                                times_dir_indep, times_dir_dep)
                            v1_phase += v1_dir_indep_interp
                        if len(times_dir_dep) != N_times_max:
                            v1_phase = interpolate_phase(v1_phase, times_dir_dep,
                                times_max)
                        phases1[i, l, k, :] = v1_phase

                    if len(times_dir_dep) != N_times_max:
                        # Set all weights to unity if interpolation was required
                        flags[i, l, k, :] = 1.0
                    else:
                        flags[i, l, k, :] = sf_dir_dep.getValues(weight=True,
                            retAxesVals=False)

                    if np.all(phases0[i, l, k, :] == 0.0):
                        # Check for flagged stations
                        flags[i, l, k, :] = 0.0
        pbar.update(i)
    pbar.finish()
    logging.info('Collection complete')
    for i, source in enumerate(source_names):
        nbands = len(find(m[i, :]))
        logging.info('  Source {0} has solutions in {1} bands'.format(source, nbands))

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

    Keyword arguments:
    phase1 -- array of phase solutions to average/interpolate
    time1 -- array of times for phase1
    time2 -- array of times to average/interpolate to
    interpMethod -- interpolation method (see scipy.interpolate.griddata)
    """
    import numpy as np
    from scipy.interpolate import interp1d

    phase1 = unwrap(phase1)

    if len(time2) < len(time1):
        # Average
        tratio = float(len(time1)) / float(len(time2))
        nslots = np.ceil(tratio)
        valsNew = average_phase(phase1, nslots)
    else:
        # Interpolate
        f = interp1d(time1, phase1, kind=interpMethod, fill_value=0.0,
            bounds_error=False)
        valsNew = f(time2)

    return valsNew


def average_phase(phase, nslots, times=None):
    """Averages input phases over nslots time slots

    Keyword arguments:
    phase -- array of phase solutions
    nslots -- average over nslots
    times -- array of times to average
    """
    import numpy as np

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
    init_sols_per_pair=False, propagate=False, nband_min=2):
    """Fits TEC values to phase solutions per source pair

    Returns TEC solutions as array of shape (N_sources, N_times, N_stations)

    Keyword arguments:
    phases -- array of phase solutions
    flags -- phase solution flags (0 = use, 1 = flagged)
    mask -- mask for sources and frequencies (0 = ignore, 1 = use)
    freqs -- array of frequencies
    init_sols -- solutions to use to initialize the fits
    init_sols_per_pair -- init_sols are per source pair (not per source)
    propagate -- propagate solutions from previous solution
    nband_min -- min number of bands for a source to be used
    """
    from pylab import pinv, newaxis, find
    import numpy as np
    from lofar.expion import baselinefitting
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar

    sols_list = []
    eq_list = []
    vars_list = []
    var_eq_list = []

    N_sources = phases.shape[0]
    N_stations = phases.shape[1]
    N_times = phases.shape[3]
    N_pairs = 0
    for i in xrange(N_sources):
        for j in xrange(i):
            subband_selection = find(mask[i, :] * mask[j, :])
            if len(subband_selection) < nband_min:
                continue
            N_pairs += 1

    if init_sols is None and not init_sols_per_pair:
        init_sols = np.zeros((N_sources, N_times, N_stations), dtype = np.float)
    elif init_sols is None and init_sols_per_pair:
        init_sols = np.zeros((N_pairs, N_times, N_stations), dtype = np.float)
    source_pairs = []

    k = 0
    ipbar = 0
    logging.info('Fitting TEC values...')
    pbar = progressbar.ProgressBar(maxval=N_pairs*N_times).start()
    for i in xrange(N_sources):
        for j in xrange(i):
            subband_selection = find(mask[i, :] * mask[j, :])
            if len(subband_selection) < nband_min:
                continue
            source_pairs.append((i, j))
            p = phases[i, :, subband_selection, :] - phases[j, :, subband_selection, :]
            A = np.zeros((len(subband_selection), 1))
            A[:, 0] = 8.44797245e9/freqs[subband_selection]

            flags_source_pair = flags[i, :, subband_selection, :] * flags[j, :,
                subband_selection, :]
            constant_parms = np.zeros((1, N_stations), dtype = np.bool)
            sols = np.zeros((N_times, N_stations), dtype = np.float)

            if init_sols_per_pair:
                p_0 = init_sols[k, 0, :][newaxis, :]
            else:
                p_0 = (init_sols[i, 0, :] - init_sols[j, 0, :])[newaxis, :]

            for t_idx in xrange(N_times):
                x = p[:, :, t_idx].copy()
                f = flags_source_pair[:, :, t_idx].copy()
                if not propagate:
                    if init_sols_per_pair:
                        p_0 = init_sols[k, t_idx, :][newaxis, :]
                    else:
                        p_0 = (init_sols[i, t_idx, :] - init_sols[j, 0, :])[newaxis, :]
                sol = baselinefitting.fit(x, A, p_0, f, constant_parms)
                if propagate:
                    p_0 = sol.copy()
                sols[t_idx, :] = sol[0, :]
                ipbar += 1
                pbar.update(ipbar)
            sols = sols[:, :] - np.mean(sols[:, :], axis=1)[:, newaxis]

            weight = len(subband_selection)
            sols_list.append(sols*weight)
            eq = np.zeros(N_sources)
            eq[i] = weight
            eq[j] = -weight
            eq_list.append(eq)
            k += 1
    pbar.finish()

    sols = np.array(sols_list)
    B = np.array(eq_list)
    source_selection = find(np.sum(abs(B), axis=0))
    N_sources = len(source_selection)
    if N_sources == 0:
        logging.error('All sources have fewer than the required minimum number of bands.')
        return None, None
    pinvB = pinv(B)
    r = np.dot(pinvB, sols.transpose([1, 0, 2]))
    r = r[source_selection, :, :]

    return r, source_selection


def add_stations(station_selection, phases0, phases1, flags, mask,
    station_names, station_positions, source_names, source_selection,
    times, freqs, r, nband_min=2, soln_type='phase', nstations_max=None,
    excluded_stations=None, t_step=5, tec_step1=5, tec_step2=21,
    search_full_tec_range=True):
    """
    Adds stations to TEC fitting using an iterative initial-guess search to
    ensure the global min is found

    Keyword arguments:
    station_selection -- indices of stations to use in fitting
    phases0 -- XX phase solutions
    phases1 -- YY phase solutions
    flags -- phase solution flags (0 = use, 1 = flagged)
    mask -- mask for sources and frequencies (0 = ignore, 1 = use)
    station_names -- array of station names
    source_names -- array of source names
    source_selection -- indices of sources to use in fitting
    times -- array of times
    freqs -- array of frequencies
    r -- array of TEC solutions returned by fit_tec_per_source_pair()
    nband_min -- min number of bands for a source to be used
    soln_type -- type of phase solution: 'phase' or 'scalarphase'
    nstations_max -- max number of stations to use
    excluded_stations -- stations to exclude
    t_step -- try full TEC range every t_step number of solution times
    tec_step1 -- number of steps in TEC subrange (+/- last TEC fit value)
    tec_step2 -- number of steps in full TEC range (-0.1 -- 0.1)
    search_full_tec_range -- always search the full TEC range (-0.1 -- 0.1)
    """
    from pylab import pinv, newaxis, find, amin
    import numpy as np
    from lofar.expion import baselinefitting
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar

    N_sources_selected = len(source_selection)
    N_stations_selected = len(station_selection)
    N_piercepoints = N_sources_selected * N_stations_selected
    N_times = len(times)
    N_stations = len(station_names)
    N_sources = len(source_names)
    N_pairs = 0
    for ii, i in enumerate(source_selection):
        for jj, j in enumerate(source_selection):
            if j == i:
                break
            subband_selection = find(mask[i,:] * mask[j,:])
            if len(subband_selection) < nband_min:
                continue
            N_pairs += 1

    D = np.resize(station_positions, (N_stations, N_stations, 3))
    D = np.transpose(D, (1, 0, 2)) - D
    D = np.sqrt(np.sum(D**2, axis=2))

    station_selection1 = station_selection
    stations_to_add = np.array([i for i in xrange(len(station_names))
        if i not in station_selection1 and station_names[i] not in
        excluded_stations])
    if len(stations_to_add) == 0:
        return station_selection1, r

    # Check if desired number of stations is already reached
    if nstations_max is not None:
        if len(station_selection1) >= nstations_max:
            return station_selection1, r

    logging.info("Using fitting with iterative search for remaining stations "
        "(up to {0} stations in total)".format(nstations_max))
    q = r
    while len(stations_to_add)>0:
        D1 = D[stations_to_add[:,newaxis], station_selection1[newaxis,:]]

        minimum_distance = amin(D1, axis=1)
        station_to_add = stations_to_add[np.argmin(minimum_distance)]
        station_selection1 = np.append(station_selection1, station_to_add)
        N_stations_selected1 = len(station_selection1)

        # Remove station from list
        stations_to_add = stations_to_add[stations_to_add != station_to_add]

        sols_list = []
        eq_list = []
        min_e_list = []

        ipbar = 0
        logging.info('Fitting TEC values with {0} included...'.format(
            station_names[station_to_add]))
        pbar = progressbar.ProgressBar(maxval=N_pairs*N_times).start()
        for ii, i in enumerate(source_selection):
            for jj, j in enumerate(source_selection):
                if j == i:
                    break
                subband_selection = find(mask[i,:] * mask[j,:])
                if len(subband_selection) < nband_min:
                    continue
                logging.debug('Adding {0} for source pair: {1}-{2}'.format(
                    station_names[station_to_add], i, j))
                p0 = phases0[i, station_selection1[:,newaxis],
                    subband_selection[newaxis,:], :] - phases0[j,
                    station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                p0 = p0 - np.mean(p0, axis=0)[newaxis,:,:]
                if soln_type != 'scalarphase':
                    p1 = phases1[i, station_selection1[:,newaxis],
                        subband_selection[newaxis,:], :] - phases1[j,
                        station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                    p1 = p1 - np.mean(p1, axis=0)[newaxis, :, :]
                A = np.zeros((len(subband_selection), 1))
                A[:, 0] = 8.44797245e9 / freqs[subband_selection]
                sd = (0.1 * 30e6) / freqs[subband_selection] # standard deviation of phase solutions as function of frequency (rad)

                flags_source_pair = flags[i, station_selection1[:,newaxis],
                    subband_selection[newaxis,:], :] * flags[j,
                    station_selection1[:,newaxis], subband_selection[newaxis,:], :]
                constant_parms = np.zeros((1, N_stations_selected1), dtype = np.bool)
                sols = np.zeros((N_times, N_stations_selected1), dtype = np.float)
                p_0_best = None
                for t_idx in xrange(N_times):
                    if np.mod(t_idx, t_step) == 0:
                        min_e = np.Inf
                        if p_0_best is not None and not search_full_tec_range:
                            min_tec = p_0_best[0, -1] - 0.02
                            max_tec = p_0_best[0, -1] + 0.02
                            nsteps = tec_step1
                        else:
                            min_tec = -0.1
                            max_tec = 0.1
                            nsteps = tec_step2
                        logging.debug('  Trying initial guesses between {0} and '
                            '{1} TECU'.format(min_tec, max_tec))
                        for offset in np.linspace(min_tec, max_tec, nsteps):
                            p_0 = np.zeros((1, N_stations_selected1), np.double)
                            p_0[0, :N_stations_selected1-1] = (q[ii, t_idx, :] -
                                q[jj, t_idx, :])[newaxis,:]
                            p_0[0, -1] = offset

                            x = p0[:,:,t_idx].copy()
                            f = flags_source_pair[:, :, t_idx].copy()
                            sol0 = baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                            sol0 -= np.mean(sol0)
                            residual = np.mod(np.dot(A, sol0) - x.T + np.pi,
                                2 * np.pi) - np.pi
                            e = np.var(residual[f.T==0])

                            if soln_type != 'scalarphase':
                                x = p1[:,:,t_idx].copy()
                                f = flags_source_pair[:, :, t_idx].copy()
                                sol1 = baselinefitting.fit(x.T, A, p_0, f,
                                    constant_parms)
                                sol1 -= np.mean(sol1)
                                residual = np.mod(np.dot(A, sol1) - x.T + np.pi,
                                    2 * np.pi) - np.pi
                                residual = residual[f.T==0]
                                e += np.var(residual)
                            else:
                                sol1 = sol0

                            if e < min_e:
                                logging.debug('  Found new min variance of {0} '
                                    'with initial guess of {1} TECU'.format(e,
                                    p_0[0, -1]))
                                min_e = e
                                p_0_best = p_0
                                sols[t_idx, :] = (sol0[0, :] + sol1[0, :])/2
                    else:
                        # Use previous init
                        x = p0[:, :, t_idx].copy()
                        f = flags_source_pair[:, :, t_idx].copy()
                        sol0 = baselinefitting.fit(x.T, A, p_0_best, f,
                            constant_parms)
                        sol0 -= np.mean(sol0)

                        if soln_type != 'scalarphase':
                            x = p1[:, :, t_idx].copy()
                            f = flags_source_pair[:, :, t_idx].copy()
                            sol1 = baselinefitting.fit(x.T, A, p_0_best, f,
                                constant_parms)
                            sol1 -= np.mean(sol1)
                        else:
                            sol1 = sol0
                        sols[t_idx, :] = (sol0[0, :] + sol1[0, :])/2

                    ipbar += 1
                    pbar.update(ipbar)

                ### Remove outliers
                logging.debug('  Searching for outliers...')
                for kk in xrange(10):
                    s = sols[:, -1].copy()
                    selection = np.zeros(len(s), np.bool)
                    for t_idx in xrange(len(s)):
                        start_idx = np.max([t_idx-10, 0])
                        end_idx = np.min([t_idx+10, len(s)])
                        selection[t_idx] = np.sum(abs(s[start_idx:end_idx] -
                            s[t_idx]) < 0.02) > (end_idx - start_idx - 8)
                    outliers = find(np.logical_not(selection))
                    if len(outliers) == 0:
                        break
                    for t_idx in outliers:
                        try:
                            idx0 = find(selection[:t_idx])[-1]
                        except IndexError:
                            idx0 = -1
                        try:
                            idx1 = find(selection[t_idx+1:])[0] + t_idx + 1
                        except IndexError:
                            idx1 = -1
                        if idx0 == -1:
                            s[t_idx] = s[idx1]
                        elif idx1 == -1:
                            s[t_idx] = s[idx0]
                        else:
                            s[t_idx] = (s[idx0] * (idx1-t_idx) + s[idx1] *
                                (t_idx-idx0)) / (idx1-idx0)

                        p_0 = np.zeros((1, N_stations_selected1), np.double)
                        p_0[0,:] = sols[t_idx,:]
                        p_0[0,-1] = s[t_idx]

                        x = p0[:,:,t_idx].copy()
                        f = flags_source_pair[:,:,t_idx].copy()
                        sol0 = baselinefitting.fit(x.T, A, p_0, f, constant_parms)
                        sol0 -= np.mean(sol0)

                        if soln_type != 'scalarphase':
                            x = p1[:,:,t_idx].copy()
                            sol1 = baselinefitting.fit(x.T, A, p_0, f,
                                constant_parms)
                            sol1 -= np.mean(sol1)
                        else:
                            sol1 = sol0
                        sols[t_idx, :] = (sol0[0,:] + sol1[0,:])/2

                weight = 1.0
                sols_list.append(weight*sols)
                min_e_list.append(min_e)
                eq = np.zeros(N_sources)
                eq[ii] = weight
                eq[jj] = -weight
                eq_list.append(eq)

        sols = np.array(sols_list)
        B = np.array(eq_list)
        pinvB = pinv(B)

        q = np.dot(pinvB, sols.transpose([1,0,2]))
        pbar.finish()

        if nstations_max is not None:
            if N_stations_selected1 == nstations_max:
                break

    return station_selection1, q


def run( step, parset, H ):
    """
    Fit phase solutions to obtain TEC values per direction.

    Phase solutions are assumed to be stored in solsets of the H5parm file, one
    solset per band per field. Only phase- or scalarphase-type solution tables
    are used. If direction-independent solutions are found (in addition to the
    direction-dependent ones), they are added, after averaging, to the
    corresponding direction-dependent ones. Phase solutions are automatically
    grouped by field and by band.

    The derived TEC values are stored in the specified output soltab of type
    'tec', with one TEC value per station per direction per solution interval.
    The TEC values are derived using the ``lofar.expion.baselinefitting.fit()``
    function to fit a TEC value to the phases. The term that is minimized
    includes all baselines, so there is no preferred reference station, and the
    residual is computed as the complex number 1.0 - exp(1i phasedifference),
    which is zero when the phase difference is a multiple of 2pi.

    The TEC solution table may be used to derive TEC screens using the
    TECSCREEN operation.
    """
    from losoto.h5parm import solFetcher, solWriter
    import numpy as np
    # Switch to the Agg backend to prevent problems with pylab imports when
    # DISPLAY env. variable is not set
    import os
    from pylab import find
    import re
    from .tecscreen import calculate_piercepoints
    from numpy.linalg import norm

    solsets = getParSolsets( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )
    nband_min = int(parset.getString('.'.join(["LoSoTo.Steps", step, "MinBands"]), '8' ))
    niter = int(parset.getString('.'.join(["LoSoTo.Steps", step, "NumIter"]), '1' ))
    dist_cut_m = np.float(parset.getString('.'.join(["LoSoTo.Steps", step, "DistCut"]), '2e3' ))
    nstations_max = int(parset.getString('.'.join(["LoSoTo.Steps", step, "MaxStations"]), '100' ))
    outSolset = parset.getString('.'.join(["LoSoTo.Steps", step, "OutSoltab"]), '' ).split('/')[0]
    outSoltab = parset.getString('.'.join(["LoSoTo.Steps", step, "OutSoltab"]), '' ).split('/')[1]

    # Collect solutions, etc. into arrays for fitting.
    (phases0, phases1, flags, mask, station_names, station_positions,
        source_names, source_positions, freqs, times, pointing,
        soln_type) = collect_solutions(H, dirs=dirs, solsets=solsets)
    if phases0 is None:
        return 1

    # Build list of stations to include
    included_stations = []
    if ants is not None:
        for ant in ants:
            included_stations += [s for s in station_names if re.search(ant, s)]
    else:
        included_stations = station_names
    excluded_stations = [s for s in station_names if s not in included_stations]

    # Select stations to use for first pass
    if len(excluded_stations) > 0:
        logging.info("Excluding stations: {0}".format(np.sort(excluded_stations)))
    mean_position = np.array([np.median(station_positions[:, 0]),
        np.median(station_positions[:, 1]), np.median(station_positions[:, 2])])
    dist = np.sqrt(np.sum((station_positions - mean_position)**2, axis=1))
    dist_sort_ind = np.argsort(dist)
    station_selection1 = find(dist < dist_cut_m)
    station_selection = np.array([i for i in dist_sort_ind
        if i in station_selection1 and station_names[i] not in excluded_stations])
    if len(station_selection) > nstations_max:
        station_selection = station_selection[:nstations_max]
    logging.info("Using normal fitting (no iterative search) for {0} stations "
        "within {1} km of the core:\n{2}".format(len(station_selection),
        dist_cut_m/1000.0, station_names[station_selection]))

    station_selection_orig = station_selection
    for iter in xrange(niter):
        # Loop over groups of nearby pierce points to identify bad stations and
        # remove them from the station_selection.
        nsig = 2.5 # number of sigma for cut
        radius = 2.0 # projected radius in km within which to compare
        if iter > 0:
            logging.info("Identifying bad stations from outlier TEC fits...")
            logging.info("Finding nearby piercepoints (assuming typical screen "
                "height of 200 km...")

            # For each source, find all the pierce points within give (projected)
            # distance from the median pierce point x,y location. Assume a typical
            # screen height of 200 km.
            pp, airmass = calculate_piercepoints(station_positions[station_selection],
                source_positions[source_selection], times, height = 200e3)
            pp = pp[0, :, :] # use first time [times, stations, dimension]
            x, y, z = station_positions[station_selection][0,:]
            east = np.array([-y, x, 0])
            east = east / norm(east)
            north = np.array([ -x, -y, (x*x + y*y)/z])
            north = north / norm(north)
            up = np.array([x ,y, z])
            up = up / norm(up)
            T = np.concatenate([east[:, np.newaxis], north[:, np.newaxis]], axis=1)
            pp1 = np.dot(pp, T).reshape((len(source_names), len(station_selection), 2))
            for i in xrange(len(source_names)):
                x_median = np.median(pp1[i, :, 0]) / 1000.0
                y_median = np.median(pp1[i, :, 1]) / 1000.0
                dist = np.sqrt( (pp1[i, :, 0] / 1000.0 - x_median)**2 +
                    (pp1[i, :, 1] / 1000.0 - y_median)**2 )
                within_radius = np.where(dist <= radius)[0]
                outside_radius = np.where(dist > radius)
                if len(within_radius) < 10:
                    logging.info("Insufficient number of closely-spaced pierce "
                        "points for bad-station detection. Skipping...")
                    abort_iter = True
                    break
                else:
                    abort_iter = False
                r_median = np.median(r[i, :, within_radius], axis=1)
                r_tot_meddiff = np.zeros(len(station_selection[within_radius]),
                    dtype=float)
                for j in xrange(len(station_selection[within_radius])):
                    r_tot_meddiff[j] = np.sum(np.abs(r[i, :, within_radius[j]] - r_median[j]))
            if abort_iter:
                break
            good_stations = np.where(r_tot_meddiff < nsig * np.median(r_tot_meddiff))
            station_selection = np.append(station_selection[within_radius[good_stations]],
                station_selection[outside_radius])

            new_excluded_stations = [station_names[s] for s in
                station_selection_orig if (s not in station_selection and
                s not in excluded_stations)]
            if len(new_excluded_stations) > 0:
                logging.info('Excluding stations due to TEC solutions that differ '
                    'significantly from mean: {0}'.format(np.sort(new_excluded_stations)))
                logging.info('Updating fit...')
                excluded_stations += new_excluded_stations
                nstations_max -= len(new_excluded_stations)
            else:
                logging.info('No bad stations found.')
                break

        # Fit a TEC value to the phase solutions per source pair.
        # No iterative search for the global minimum is done
        if soln_type == 'scalarphase':
            r, source_selection = fit_tec_per_source_pair(
                phases0[:, station_selection, :, :],
                flags[:, station_selection, :, :],
                mask, freqs, propagate=True, nband_min=nband_min)
            if r is None:
                return 1
        else:
            r0, source_selection = fit_tec_per_source_pair(
                phases0[:, station_selection, :, :],
                flags[:, station_selection, :, :],
                mask, freqs, propagate=True, nband_min=nband_min)
            r1, source_selection = fit_tec_per_source_pair(
                phases1[:, station_selection, :, :],
                flags[:, station_selection, :, :],
                mask, freqs, propagate=True, nband_min=nband_min)
            if r0 is None or r1 is None:
                return 1

            # take the mean of the two polarizations
            r = (r0 + r1) / 2

    # Add stations by searching iteratively for global minimum in solution space
    station_selection, r = add_stations(station_selection, phases0,
        phases1, flags, mask, station_names, station_positions, source_names,
        source_selection, times, freqs, r, nband_min=nband_min,
        soln_type=soln_type, nstations_max=nstations_max, excluded_stations=
        excluded_stations, search_full_tec_range=False)

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

    tf_st = H.makeSoltab(solset._v_name, 'tec', outSoltab,
        axesNames=['dir', 'time', 'ant'], axesVals=[dirs_out, times, ants_out],
        vals=r[source_selection, :, :],
        weights=np.ones_like(r[source_selection, :, :]))

    # Add history
    sw = solWriter(tf_st)
    sw.addHistory('CREATE (by TECFIT operation)')

    return 0


