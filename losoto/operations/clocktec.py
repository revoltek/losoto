#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Clock/TEC separation operation for LoSoTo, based on the lofar_clocktec package
# (pip install lofar-clocktec).

from losoto._logging import logger as logging

logging.debug('Loading CLOCKTEC module.')

# Polarizations that carry the clock/TEC signal (cross-hands are skipped).
_PARALLEL_POLS = ['XX', 'YY', 'RR', 'LL', 'I']


def _run_parser(soltab, parser, step):
    clocksoltabOut = parser.getstr( step, 'clocksoltabOut', 'clock000' )
    tecsoltabOut = parser.getstr( step, 'tecsoltabOut', 'tec000' )
    offsetsoltabOut = parser.getstr( step, 'offsetsoltabOut', 'phase_offset000' )
    tec3rdsoltabOut = parser.getstr( step, 'tec3rdsoltabOut', 'tec3rd000' )
    refAnt = parser.getstr( step, 'refAnt', '' )
    numParams = parser.getint( step, 'numParams', 0 )
    combinePol = parser.getbool( step, 'combinePol', False )
    circular = parser.getbool( step, 'circular', False )
    weightMode = parser.getstr( step, 'weightMode', 'kappa' )
    initMode = parser.getstr( step, 'initMode', 'robust' )
    initWithPrevious = parser.getbool( step, 'initWithPrevious', True )
    removeSteps = parser.getbool( step, 'removeSteps', False )
    spatialFit = parser.getbool( step, 'spatialFit', False )
    allowPiWraps = parser.getbool( step, 'allowPiWraps', False )
    smoothSolutions = parser.getbool( step, 'smoothSolutions', False )
    lofar2Mode = parser.getbool( step, 'lofar2Mode', True )
    finalPh0Fit = parser.getint( step, 'finalPh0Fit', 0 )
    logFile = parser.getstr( step, 'logFile', '' )

    parser.checkSpelling( step, soltab, ['tecsoltabOut', 'clocksoltabOut', 'offsetsoltabOut', 'tec3rdsoltabOut',
                                         'refAnt', 'numParams', 'combinePol', 'circular',
                                         'weightMode', 'initMode', 'initWithPrevious', 'removeSteps', 'spatialFit',
                                         'allowPiWraps', 'smoothSolutions', 'lofar2Mode', 'finalPh0Fit', 'logFile'] )

    return run( soltab, tecsoltabOut=tecsoltabOut, clocksoltabOut=clocksoltabOut, offsetsoltabOut=offsetsoltabOut,
                tec3rdsoltabOut=tec3rdsoltabOut, refAnt=refAnt, numParams=numParams,
                combinePol=combinePol, circular=circular, weightMode=weightMode, initMode=initMode,
                initWithPrevious=initWithPrevious, removeSteps=removeSteps, spatialFit=spatialFit,
                allowPiWraps=allowPiWraps, smoothSolutions=smoothSolutions, lofar2Mode=lofar2Mode,
                finalPh0Fit=finalPh0Fit, logFile=logFile )


def _select_ref_station(stations, refAnt):
    """Return the index of the reference station."""
    stations = list(stations)
    if refAnt:
        if refAnt not in stations:
            raise ValueError('Reference antenna %s not found in the selected antennas.' % refAnt)
        return stations.index(refAnt)
    # default: CS002 (centre of the superterp), otherwise the superterp station, otherwise the first one
    for prefix in ['CS002', 'ST001']:
        for i, st in enumerate(stations):
            if st.startswith(prefix):
                return i
    return 0


def _fit_one(phases, times, freqs, stations, fit_kwargs):
    """
    Run lofar_clocktec on a single (time, freq, ant) phase cube.
    Flagged samples must be NaN in phases.

    Returns clock [ns], tec [TECU], phase0 [rad], tec3rd, all with shape (time, ant).
    """
    import numpy as np
    import astropy.units as u
    from astropy.time import Time
    from lofar_clocktec import CTInput, fit_clock_tec

    ct_input = CTInput(
        phases=phases * u.rad,
        freqs=np.asarray(freqs, dtype=float) * u.Hz,
        flags=np.isnan(phases),  # informative only, the fitter uses NaNs
        stations=[str(s) for s in stations],
        times=Time(np.asarray(times, dtype=float) / 86400., format='mjd'),  # LoSoTo times are MJD in seconds
    )
    result = fit_clock_tec(ct_input, **fit_kwargs)

    # lofar_clocktec returns arrays of shape (ant, time): transpose to (time, ant)
    def _tr(arr):
        if arr is None:
            return np.full((len(times), len(stations)), np.nan)
        return np.asarray(arr, dtype=float).T

    return _tr(result.clock), _tr(result.tec), _tr(result.phase0), _tr(result.tec3rd)


def run( soltab, tecsoltabOut='tec000', clocksoltabOut='clock000', offsetsoltabOut='phase_offset000',
         tec3rdsoltabOut='tec3rd000', refAnt='', numParams=0, combinePol=False, circular=False,
         weightMode='kappa', initMode='robust', initWithPrevious=True, removeSteps=False, spatialFit=False,
         allowPiWraps=False, smoothSolutions=False, lofar2Mode=True, finalPh0Fit=0, logFile='' ):
    """
    Separate phase solutions into Clock and TEC using the lofar_clocktec package.
    The Clock, TEC, phase offset (and optionally 3rd order TEC) values are stored in the
    specified output soltabs with type 'clock', 'tec', 'phase' and 'tec3rd'.

    Parameters
    ----------
    tecsoltabOut : str, optional
        Name of the output TEC soltab, by default 'tec000'.

    clocksoltabOut : str, optional
        Name of the output clock soltab (in seconds), by default 'clock000'.

    offsetsoltabOut : str, optional
        Name of the output phase offset soltab (in radians), by default 'phase_offset000'.
        Contrary to the old implementation it has a time axis.

    tec3rdsoltabOut : str, optional
        Name of the output 3rd order TEC soltab, by default 'tec3rd000'. Only written if numParams is 4.

    refAnt : str, optional
        Reference antenna. All phases are referenced to it before fitting, so its solutions are zero.
        By default the first CS002 station (or ST001, or the first antenna).

    numParams : int, optional
        Number of fitted parameters: 2 (clock + TEC), 3 (+ phase offset), 4 (+ 3rd order TEC,
        useful < 40 MHz).
        By default 0: automatic choice (2 for HBA or narrow-band data, 3 otherwise).

    combinePol : bool, optional
        Find a combined polarization solution, by default False.

    circular : bool, optional
        Assume circular polarization with FR not removed. Only relevant with combinePol=True:
        the RR and LL phases are summed (cancelling Faraday rotation) and the results halved.
        By default False.

    weightMode : str, optional
        Per-channel weighting: 'kappa' (Von Mises concentration, robust) or 'complex'. By default 'kappa'.

    initMode : str, optional
        Initialisation strategy: 'robust' (brute-force phasor search) or 'linear' (faster). By default 'robust'.

    initWithPrevious : bool, optional
        Initialise each timestep with the previous solution (much faster). By default True.

    removeSteps : bool, optional
        Correct clock/TEC ambiguity jumps between timesteps (recommended for HBA). By default False.

    spatialFit : bool, optional
        Fit a spatial polynomial to the TEC to fix integer ambiguities and derive phase offsets.
        Recommended when phase offsets are needed. By default False.

    allowPiWraps : bool, optional
        Allow pi phase jumps (e.g. international stations with imperfect Faraday correction). By default False.

    smoothSolutions : bool, optional
        Smooth clock and TEC solutions in time after fitting. By default False.

    lofar2Mode : bool, optional
        Restrict the delay search range to LOFAR2 timing accuracy (<20 ns). Set to False for LOFAR1 data
        where remote station clocks can be much larger. By default True.

    finalPh0Fit : int, optional
        Experimental constant phase offset fit: 0 (off), 1 (constant ph0), 2 (constant ph0 and clock).
        By default 0.

    logFile : str, optional
        Write the detailed lofar_clocktec log to this file. By default '' (use the Python logger).

    Notes
    -----
    The phase offsets are stored with the same sign convention as clock and TEC (part of the
    model phase), so they can be subtracted directly with the residuals operation.
    """
    import numpy as np

    try:
        import lofar_clocktec  # noqa: F401
    except ImportError:
        logging.error('The CLOCKTEC operation requires the lofar_clocktec package (pip install lofar-clocktec).')
        return 1

    logging.info("Clock/TEC separation on soltab: "+soltab.name)

    # some checks
    solType = soltab.getType()
    if solType != 'phase':
        logging.warning("Soltab type of "+soltab.name+" is: "+solType+" should be phase. Ignoring.")
        return 1
    if weightMode not in ['kappa', 'complex']:
        logging.error('weightMode must be "kappa" or "complex".')
        return 1
    if initMode not in ['robust', 'linear']:
        logging.error('initMode must be "robust" or "linear".')
        return 1
    if numParams not in [0, 2, 3, 4]:
        logging.error('numParams must be 0 (auto), 2, 3 or 4.')
        return 1
    do3rd = (numParams == 4)

    axesNames = soltab.getAxesNames()
    for ax in ['ant', 'freq', 'time']:
        if ax not in axesNames:
            logging.error('Clock/TEC separation needs a soltab with a "%s" axis.' % ax)
            return 1
    hasPol = 'pol' in axesNames

    stations = soltab.getAxisValues('ant')
    freqs = soltab.getAxisValues('freq')
    times = soltab.getAxisValues('time')
    if len(stations) < 2:
        logging.error('Clock/TEC separation needs at least 2 antennas selected.')
        return 1
    if len(freqs) < 10:
        logging.error('Clock/TEC separation needs at least 10 frequency channels, preferably distributed over a wide range')
        return 1

    try:
        refIdx = _select_ref_station(stations, refAnt)
    except ValueError as e:
        logging.error(str(e))
        return 1
    logging.info('Using %s as reference antenna.' % stations[refIdx])

    # polarizations to fit (parallel hands only)
    if hasPol:
        pols = list(soltab.getAxisValues('pol'))
        polIdx = [i for i, p in enumerate(pols) if p in _PARALLEL_POLS]
        if len(polIdx) == 0:
            polIdx = list(range(len(pols)))
        fitPols = [pols[i] for i in polIdx]
        if combinePol and len(polIdx) != 2:
            logging.error('combinePol needs exactly two parallel-hand polarizations, found: %s' % fitPols)
            return 1
        if len(polIdx) < len(pols):
            logging.info('Fitting polarizations: %s' % fitPols)
    circularCombine = combinePol and circular

    fit_kwargs = dict(weight_mode=weightMode, initialize_mode=initMode, init_with_previous=initWithPrevious,
                      num_params_init=(numParams if numParams else None), final_ph0_fit=finalPh0Fit,
                      remove_steps=removeSteps, spatial_fit=spatialFit, allow_pi_wraps=allowPiWraps,
                      ref_station=refIdx, smooth_solutions=smoothSolutions, lofar2_mode=lofar2Mode,
                      log_file=(logFile if logFile else None))

    # output axes: all input axes but freq (and pol if combined), time and ant first
    otherAxes = [ax for ax in axesNames if ax not in ['time', 'ant', 'freq', 'pol']]
    outAxes = ['time', 'ant'] + (['pol'] if hasPol and not combinePol else []) + otherAxes
    outVals = {'time': times, 'ant': stations}
    if hasPol and not combinePol:
        outVals['pol'] = fitPols
    for ax in otherAxes:
        outVals[ax] = soltab.getAxisValues(ax)
    outShape = [len(outVals[ax]) for ax in outAxes]
    clockOut = np.zeros(outShape, dtype=float)
    tecOut = np.zeros(outShape, dtype=float)
    offsetOut = np.zeros(outShape, dtype=float)
    tec3rdOut = np.zeros(outShape, dtype=float)
    weightsOut = np.zeros(outShape, dtype=float)

    returnAxes = ['time', 'freq', 'ant'] + (['pol'] if hasPol else [])
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=returnAxes, weight=True):

        # reorder to (time, freq, ant[, pol]) and flag with NaNs
        iterAxes = [ax for ax in axesNames if ax in returnAxes]
        order = [iterAxes.index(ax) for ax in returnAxes]
        vals = np.array(vals, dtype=float).transpose(order)
        weights = np.array(weights, dtype=float).transpose(order)
        vals[weights == 0] = np.nan
        if not hasPol:
            vals = vals[..., np.newaxis]
        else:
            vals = vals[..., polIdx]

        if combinePol:
            if circularCombine:
                # RR + LL cancels Faraday rotation; the doubled clock/TEC is halved afterwards
                vals = np.sum(vals, axis=-1, keepdims=True)
            else:
                phasors = np.exp(1j * vals)
                allFlagged = np.all(np.isnan(vals), axis=-1, keepdims=True)
                vals = np.angle(np.nansum(phasors, axis=-1, keepdims=True))
                vals[allFlagged] = np.nan

        # position of this iteration in the output arrays
        otherIdx = tuple(list(outVals[ax]).index(coord[ax]) for ax in otherAxes)

        for ipol in range(vals.shape[-1]):
            if hasPol and not combinePol:
                logging.info('Fitting polarization %s' % fitPols[ipol])
            clock, tec, offset, tec3rd = _fit_one(vals[..., ipol], times, freqs, stations, fit_kwargs)

            if circularCombine:
                clock, tec, offset, tec3rd = clock / 2., tec / 2., offset / 2., tec3rd / 2.

            valid = np.isfinite(clock) & np.isfinite(tec)
            valid[:, refIdx] = True  # reference station is zero by definition
            if hasPol and not combinePol:
                idx = (slice(None), slice(None), ipol) + otherIdx
            else:
                idx = (slice(None), slice(None)) + otherIdx
            clockOut[idx] = np.where(valid, clock, 0.)
            tecOut[idx] = np.where(valid, tec, 0.)
            offsetOut[idx] = np.where(valid & np.isfinite(offset), offset, 0.)
            tec3rdOut[idx] = np.where(valid & np.isfinite(tec3rd), tec3rd, 0.)
            weightsOut[idx] = valid.astype(float)

            nbad = np.sum(~valid)
            if nbad > 0:
                logging.info('Flagged %i of %i clock/TEC solutions.' % (nbad, valid.size))

    solset = soltab.getSolset()
    axesVals = [outVals[ax] for ax in outAxes]
    st = solset.makeSoltab('tec', soltabName=tecsoltabOut, axesNames=outAxes, axesVals=axesVals,
                           vals=tecOut, weights=weightsOut)
    st.addHistory('CREATE (by CLOCKTEC operation, lofar_clocktec)')
    st = solset.makeSoltab('clock', soltabName=clocksoltabOut, axesNames=outAxes, axesVals=axesVals,
                           vals=clockOut*1e-9, weights=weightsOut)  # ns -> s
    st.addHistory('CREATE (by CLOCKTEC operation, lofar_clocktec)')
    st = solset.makeSoltab('phase', soltabName=offsetsoltabOut, axesNames=outAxes, axesVals=axesVals,
                           vals=offsetOut, weights=weightsOut)
    st.addHistory('CREATE (by CLOCKTEC operation, lofar_clocktec)')
    if do3rd:
        st = solset.makeSoltab('tec3rd', soltabName=tec3rdsoltabOut, axesNames=outAxes, axesVals=axesVals,
                               vals=tec3rdOut, weights=weightsOut)
        st.addHistory('CREATE (by CLOCKTEC operation, lofar_clocktec)')

    return 0
