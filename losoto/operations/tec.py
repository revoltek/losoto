#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp
import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d
from losoto.lib_unwrap import unwrap_2d
from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading TEC module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'tec000' )
    refAnt = parser.getstr( step, 'refAnt', '')
    maxResidualFlag = parser.getfloat( step, 'maxResidualFlag', 2.5 )
    maxResidualProp = parser.getfloat( step, 'maxResidualProp', 1. )
    ncpu = parser.getint( '_global', 'ncpu', 0 )

    parser.checkSpelling( step, soltab, ['soltabOut', 'refAnt', 'maxResidualFlag', 'maxResidualProp'])
    return run(soltab, soltabOut, refAnt, maxResidualFlag, maxResidualProp, ncpu)


def mod(d):
    """ wrap phases to (-pi,pi)"""
    return np.mod(d + np.pi, 2. * np.pi) - np.pi

_noiseweight = None


def tec_merit_brute(dTEC, freq, phases, weight=True):
    """
    Merit function for brute-force grid search.
    Parameters
    ----------
    dTEC: float, dTEC in TECU
    freq: np array of floats, phases in Hz
    phases: phases to fit in rad
    weight: bool, default True. Whether to weight residuals by noise or not

    Returns
    -------
    rms_phase_residual: float, rms of phase residuals
    """
    if weight:
        w = _noiseweight(freq)
        w /= np.mean(w)
        rms_phase_residual = np.std(w * (mod(-8.44797245e9 * dTEC / freq - phases)))
    else:
        rms_phase_residual = np.std(mod(-8.44797245e9 * dTEC / freq - phases))
    return rms_phase_residual


def tec_merit(dTEC, freq, phases, weight=True):
    """
    Merit function for least-square fit.
    Parameters
    ----------
    dTEC: float, dTEC in TECU
    freq: np array of floats, phases in Hz
    phases: phases to fit in rad
    weight: bool, default True. Whether to weight residuals by noise or not

    Returns
    -------
    rms_phase_residuals: array of floats, rms of phase residuals
    """
    if weight:
        w = _noiseweight(freq)
        w /= np.mean(w)
        rms_phase_residuals = w * mod(-8.44797245e9 * dTEC / freq - phases)
    else:
        rms_phase_residuals = mod(-8.44797245e9 * dTEC / freq - phases)
    return rms_phase_residuals


def fit_tec_to_phases(vals, weights, coord, refAnt, maxResidualFlag, maxResidualProp):
    """
    Fit dTEC to phases and frequencies for a range of times.
    Parameters
    ----------
    vals: numpy array of floats, phases to fit in rad. Shape (n_times, n_freq)
    weights: numpy array of floats, phases to fit in rad. Shape (n_times, n_freq)
    coord: dict of coords of current selection. Contains time, freq, ant, (optional: dir)
    refAnt: string, reference antenna
    maxResidualFlag: float, default = 2.5 Maximum residual that is not flagged. 0=don't flag.
    maxResidualProp: float, default = 1. Maximum residual that is propagated. 0=propagate all.

    Returns
    -------
    [fitdTEC, fitweights]: list of numpy array of floats, dTEC restults in TECU / weights
    """
    # Prepare output arrays
    fitdTEC = np.zeros(len(coord['time']))
    fitweights = np.ones(len(coord['time']))  # all unflagged to start
    # skip refAnt
    if coord['ant'] == refAnt:
        return fitdTEC, fitweights

    # find flagged times, either fully flagged or less than 10 freq points...
    flagged_t = np.sum(weights, axis=1)
    flagged_t = flagged_t < 10
    flag_frac = np.sum(flagged_t) / len(flagged_t)
    if flag_frac > 0.1:
        logging.info(f'''Times with less than 10 unflagged freqs: {coord['time'][flagged_t]}: percentage: {flag_frac:.1%}''')

    ranges, Ns = (-0.5, 0.5), 1000  # default range for brute grid minimization, size of grid
    freq = coord['freq'].copy()
    # Iterate times
    for t, (time, phases, w_t) in enumerate(zip(coord['time'],vals,weights)):
        w_t = w_t.astype(bool)
        if sum(w_t) < 10:
            fitdTEC[t] = 0.
            fitweights[t] = 0
            continue
        # brute force to find global minimum
        dTEC_gridsearch = scipy.optimize.brute(tec_merit_brute, ranges=(ranges,), Ns=Ns, args=(freq[w_t], phases[w_t]))[0]
        result, success = scipy.optimize.leastsq(tec_merit, dTEC_gridsearch, args=(freq[w_t], phases[w_t]))
        best_residual = tec_merit_brute(result, freq[w_t], phases[w_t])
        # logging.info(f'result {result} cost {best_residual}')
        fitdTEC[t] = result

        if maxResidualFlag == 0 or best_residual < maxResidualFlag:
            fitweights[t] = 1
            if maxResidualProp == 0 or best_residual < maxResidualProp:
                ranges = (fitdTEC[t] - 0.05, fitdTEC[t] + 0.05)
                Ns = 100
            else:
                ranges = (-0.5, 0.5)
                Ns = 1000
        else:
            # high residual, flag and reset initial guess
            if 'dir' in coord.keys():
                logging.warning('Bad solution for ant: ' + coord['ant'] + '; dir: ' + coord['dir'] + ' (time: ' + str(
                    t) + ', resdiual: ' + str(best_residual) + ').')
            else:
                logging.warning('Bad solution for ant: ' + coord['ant'] + ' (time: ' + str(t) + ', resdiual: ' + str(
                    best_residual) + ').')
            fitweights[t] = 0
            ranges = (-0.5, 0.5)
            Ns = 1000

            # Debug plot
            # doplot = False
            # if doplot and (coord['ant'] == 'RS509LBA' or coord['ant'] == 'RS210LBA') and t % 50 == 0:
            #     print("Plotting")
            #     if not 'matplotlib' in sys.modules:
            #         import matplotlib as mpl
            #
            #         mpl.rc('figure.subplot', left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.22, hspace=0.22)
            #         mpl.use("Agg")
            #     import matplotlib.pyplot as plt
            #
            #     fig = plt.figure()
            #     fig.subplots_adjust(wspace=0)
            #     ax = fig.add_subplot(111)
            #
            #     # plot rm fit
            #     plotd = lambda d, freq: -8.44797245e9 * d / freq
            #     ax.plot(freq, plotd(fitresultd[0], freq[:]), "-", color='purple')
            #     ax.plot(freq, mod(plotd(fitresultd[0], freq[:])), ":", color='purple')
            #
            #     # ax.plot(freq, vals[idx,t], '.b' )
            #     # ax.plot(freq, phaseComb + numjumps * 2*np.pi, 'x', color='purple' )
            #     ax.plot(freq, phaseComb, 'o', color='purple')
            #
            #     residual = mod(plotd(fitd[t], freq[:]) - phaseComb)
            #     ax.plot(freq, residual, '.', color='orange')
            #
            #     ax.set_xlabel('freq')
            #     ax.set_ylabel('phase')
            #     # ax.set_ylim(ymin=-np.pi, ymax=np.pi)
            #
            #     logging.warning('Save pic: ' + str(t) + '_' + coord['ant'] + '.png')
            #     plt.savefig(str(t) + '_' + coord['ant'] + '.png', bbox_inches='tight')
            #     del fig
    if 'dir' in coord.keys():
        logging.info('%s; %s: average tec: %f TECU; std tec: %f TECU' % (coord['ant'], coord['dir'], np.mean(fitdTEC), np.std(fitdTEC))) # prev. factor of 2? why?
    else:
        logging.info('%s: average tec: %f TECU; std tec: %f TECU' % (coord['ant'], np.mean(fitdTEC), np.std(fitdTEC)))
    return [fitdTEC, fitweights]


def run( soltab, soltabOut, refAnt, maxResidualFlag, maxResidualProp, ncpu ):
    """
    Bruteforce TEC extraction from phase solutions.

    Parameters
    ----------
    soltabOut : str, optional
        output table name (same solset), by deault "tec".

    refAnt : str, optional
        Reference antenna, by default the first.

    maxResidualFlag : float, optional
        Max average residual in radians before flagging datapoint, by default 2.5 If 0: no check.

    maxResidualProp : float, optional
        Max average residual in radians before stop propagating solutions, by default 1. If 0: no check.
    """
    # define weight function for fitting (global scope)
    pth = os.path.dirname(os.path.abspath(__file__))
    wfreq, sefd = np.loadtxt(pth+'/../data/SEFD_LBA_ALL.csv', delimiter=',').T
    wgt = 1 / sefd
    wgt /= wgt.max()
    global _noiseweight
    _noiseweight = interp1d(wfreq, wgt, bounds_error=False, fill_value=(0, 1))

    logging.info("Find TEC for soltab: "+soltab.name)

    # input check
    solType = soltab.getType()
    if solType != 'phase':
       logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase. Ignoring.")
       return 1

    if 'pol' in soltab.getAxesNames():
        logging.warning("Soltab with pol axis not supported for TEC extraction. Ignoring.")
        return 1
    ants = soltab.getAxisValues('ant')
    if refAnt != '' and refAnt != 'closest' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.error('Reference antenna '+refAnt+' not found. Using: '+ants[0])
        refAnt = ants[0]
    if refAnt == '': refAnt = ants[0]

    # times and ants needs to be complete or selection is much slower
    times = soltab.getAxisValues('time')

    # create new table
    axes = soltab.getAxesNames() # ['time','freq','ant'] or ['time','freq','ant','dir']
    outaxes = axes.copy()
    outaxes.remove('freq') # ['time','ant'] or ['time','ant','dir']
    results_dTEC = np.zeros(shape=tuple(soltab.getAxisLen(axisName) for axisName in outaxes))
    results_w =  np.ones(shape=tuple(soltab.getAxisLen(axisName) for axisName in outaxes))
    solset = soltab.getSolset()
    if 'tec000' in solset.getSoltabNames():
        logging.warning('Soltab tec000 exists. Overwriting...')
        solset.getSoltab('tec000').delete()
    soltabout = solset.makeSoltab(soltype = 'tec', soltabName = soltabOut, axesNames=outaxes, \
                                  axesVals=[soltab.getAxisValues(axisName) for axisName in outaxes], \
                      vals=np.zeros(shape=tuple(soltab.getAxisLen(axisName) for axisName in outaxes)), \
                      weights=np.ones(shape=tuple(soltab.getAxisLen(axisName) for axisName in outaxes)) )
    soltabout.addHistory('Created by TEC operation from %s.' % soltab.name)

    # Collect arguments for pool.map()
    args = []
    selections = []
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['time','freq'], weight=True, refAnt=refAnt):
        if len(coord['freq']) < 10:
            logging.error('Delay estimation needs at least 10 frequency channels, preferably distributed over a wide range.')
            return 1
        args.append([vals, weights, coord, refAnt, maxResidualFlag, maxResidualProp])
        selections.append(selection)

    if ncpu == 0:
        ncpu = mp.cpu_count()
    with mp.Pool(ncpu) as pool:
        logging.info('Start TEC fitting.')
        results = pool.starmap(fit_tec_to_phases, args)

    # reorder results
    for selection, result in zip(selections,results):
        selection = [selection[0],*selection[2:]]
        print(selection,results_dTEC.shape)
        results_dTEC[selection] = np.resize(result[0], results_dTEC[selection].shape)
        results_w[selection] = np.resize(result[1], results_dTEC[selection].shape)
    # write results
    soltabout.setValues( results_dTEC )
    soltabout.setValues( results_w, weight=True )

    return 0
