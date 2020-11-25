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
    # skip refAnt and ants in stationconstraint with refAnt
    if coord['ant'] == refAnt or np.all(vals == 0.):
        pass
    else:
        # find flagged times, either fully flagged or less than 10 freq points...
        flagged_t = np.sum(weights, axis=1)
        flagged_t = flagged_t < 10
        flag_frac = np.sum(flagged_t) / len(flagged_t)
        if flag_frac > 0.1:
            logging.warning(f'''Times with less than 10 unflagged freqs: {coord['time'][flagged_t]}: percentage: {flag_frac:.1%}''')

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
    # TODO: does not work for global selection, i.e. ant=[...] in parset.
    # define weight function for fitting (global scope)
    global _noiseweight
    # Freq/SEFD estimate for LBA lofar2.0, TODO: move to data dir, add lofar outer/inner
    data_lba_all = np.array([[2.850000000000000000e+07, 2.994030080626945710e+04],
                             [2.903030303030303121e+07, 2.913874315792467314e+04],
                             [2.956060606060606241e+07, 2.848637311863995637e+04],
                             [3.009090909090908989e+07, 2.793457127624304485e+04],
                             [3.062121212121212110e+07, 2.726112344516822850e+04],
                             [3.115151515151515231e+07, 2.656675093925880356e+04],
                             [3.168181818181817979e+07, 2.597196672786273120e+04],
                             [3.221212121212121099e+07, 2.571364929933567691e+04],
                             [3.274242424242424220e+07, 2.546725999266218423e+04],
                             [3.327272727272727340e+07, 2.514084102040380822e+04],
                             [3.380303030303030461e+07, 2.486516674111211978e+04],
                             [3.433333333333333582e+07, 2.467008078555054453e+04],
                             [3.486363636363635957e+07, 2.456891934402090556e+04],
                             [3.539393939393939078e+07, 2.449906041042739525e+04],
                             [3.592424242424242198e+07, 2.441152557008936856e+04],
                             [3.645454545454545319e+07, 2.435284099484204125e+04],
                             [3.698484848484848440e+07, 2.428130362751187931e+04],
                             [3.751515151515151560e+07, 2.412378187290063579e+04],
                             [3.804545454545454681e+07, 2.398172239664983135e+04],
                             [3.857575757575757802e+07, 2.395988451670273935e+04],
                             [3.910606060606060922e+07, 2.430045490772785342e+04],
                             [3.963636363636364043e+07, 2.471031800485042186e+04],
                             [4.016666666666666418e+07, 2.386110364099146318e+04],
                             [4.069696969696969539e+07, 2.353690296084890724e+04],
                             [4.122727272727272660e+07, 2.336803344857854609e+04],
                             [4.175757575757575780e+07, 2.337549122132260163e+04],
                             [4.228787878787878901e+07, 2.282559852902076818e+04],
                             [4.281818181818181276e+07, 2.288857468625183537e+04],
                             [4.334848484848484397e+07, 2.293901290797947513e+04],
                             [4.387878787878787518e+07, 2.279201399700451293e+04],
                             [4.440909090909090638e+07, 2.263575546276593013e+04],
                             [4.493939393939393759e+07, 2.234381817033494008e+04],
                             [4.546969696969696879e+07, 2.201090887195047253e+04],
                             [4.600000000000000000e+07, 2.167511043603119106e+04],
                             [4.653030303030303121e+07, 2.126023280379924108e+04],
                             [4.706060606060606241e+07, 2.094471500899021703e+04],
                             [4.759090909090909362e+07, 2.076310086716620208e+04],
                             [4.812121212121212482e+07, 2.072452147961255469e+04],
                             [4.865151515151514858e+07, 2.072524703336137827e+04],
                             [4.918181818181817979e+07, 2.024623970337267383e+04],
                             [4.971212121212121099e+07, 2.024902352195448111e+04],
                             [5.024242424242424220e+07, 2.042157956576801007e+04],
                             [5.077272727272727340e+07, 2.061233436890421581e+04],
                             [5.130303030303029716e+07, 2.052460108698517070e+04],
                             [5.183333333333332837e+07, 2.044476653198705753e+04],
                             [5.236363636363635957e+07, 2.029609696197773519e+04],
                             [5.289393939393939078e+07, 2.011980530139794064e+04],
                             [5.342424242424242198e+07, 2.000652821604120982e+04],
                             [5.395454545454545319e+07, 1.972491939080056181e+04],
                             [5.448484848484848440e+07, 1.970762288441187411e+04],
                             [5.501515151515151560e+07, 1.963280298568139551e+04],
                             [5.554545454545454681e+07, 1.947178074763434779e+04],
                             [5.607575757575757802e+07, 1.930311130654076624e+04],
                             [5.660606060606060922e+07, 1.930880388910132751e+04],
                             [5.713636363636363298e+07, 1.896554359589192609e+04],
                             [5.766666666666666418e+07, 1.897491931565966297e+04],
                             [5.819696969696969539e+07, 1.899107079252671974e+04],
                             [5.872727272727272660e+07, 1.897387182282819413e+04],
                             [5.925757575757575780e+07, 1.910503545843373649e+04],
                             [5.978787878787878156e+07, 1.940954548637071275e+04],
                             [6.031818181818181276e+07, 1.984135252611626129e+04],
                             [6.084848484848484397e+07, 2.014750239689876980e+04],
                             [6.137878787878787518e+07, 2.045096111245894645e+04],
                             [6.190909090909090638e+07, 2.102328279043795555e+04],
                             [6.243939393939393759e+07, 2.125787547095409172e+04],
                             [6.296969696969696879e+07, 2.165665026710620077e+04],
                             [6.350000000000000000e+07, 2.219586934836382716e+04],
                             [6.403030303030303121e+07, 2.251618565684030182e+04],
                             [6.456060606060606241e+07, 2.256982609767890244e+04],
                             [6.509090909090908617e+07, 2.299487046699088387e+04],
                             [6.562121212121211737e+07, 2.360399101053923732e+04],
                             [6.615151515151514858e+07, 2.452839877523585892e+04],
                             [6.668181818181817979e+07, 2.475335631368393661e+04],
                             [6.721212121212121844e+07, 2.541818455542022275e+04],
                             [6.774242424242424965e+07, 2.632339058027763167e+04],
                             [6.827272727272728086e+07, 2.722260078658507700e+04],
                             [6.880303030303029716e+07, 2.809248882867760403e+04],
                             [6.933333333333332837e+07, 2.869324085220565030e+04],
                             [6.986363636363635957e+07, 2.937062652000980597e+04],
                             [7.039393939393939078e+07, 2.975876253021111552e+04],
                             [7.092424242424242198e+07, 3.023381387615679705e+04],
                             [7.145454545454545319e+07, 3.078095987549442361e+04],
                             [7.198484848484848440e+07, 3.136289409735724985e+04],
                             [7.251515151515151560e+07, 3.194191655745387106e+04],
                             [7.304545454545454681e+07, 3.256773533922632851e+04],
                             [7.357575757575756311e+07, 3.323538538253790466e+04],
                             [7.410606060606059432e+07, 3.408304345233587082e+04],
                             [7.463636363636362553e+07, 3.532438151289376401e+04],
                             [7.516666666666665673e+07, 3.611054891041758674e+04],
                             [7.569696969696968794e+07, 3.677927807149000728e+04],
                             [7.622727272727271914e+07, 3.742928695890052768e+04],
                             [7.675757575757575035e+07, 3.831019856693004840e+04],
                             [7.728787878787878156e+07, 3.942348853846100246e+04],
                             [7.781818181818181276e+07, 3.970818633546015189e+04],
                             [7.834848484848484397e+07, 4.055220533173037984e+04],
                             [7.887878787878787518e+07, 4.174217174788207922e+04],
                             [7.940909090909090638e+07, 4.371878037581651733e+04],
                             [7.993939393939393759e+07, 4.366688385026387550e+04],
                             [8.046969696969696879e+07, 4.441932939676849492e+04],
                             [8.100000000000000000e+07, 4.513826943235301587e+04]])
    sefd_lba_all, wfreq = data_lba_all[:,1], data_lba_all[:,0]
    wgt = 1 / sefd_lba_all
    wgt /= wgt.max()
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
        logging.warning('Reference antenna '+refAnt+' not found. Using: '+ants[0])
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
        selection = tuple([axsel for axsel, ax in zip(selection,axes) if ax in outaxes]) # get rid of selection along freq axis
        results_dTEC[selection] = np.resize(result[0], results_dTEC[selection].shape)
        results_w[selection] = np.resize(result[1], results_dTEC[selection].shape)
    # write results
    soltabout.setValues( results_dTEC )
    soltabout.setValues( results_w, weight=True )

    return 0
