#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading POLALIGN module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'phasediff' )
    minFreq = parser.getfloat( step, 'minFreq', 0 )
    refAnt = parser.getstr( step, 'refAnt', '' )

    parser.checkSpelling( step, soltab, ['soltabOut', 'minFreq', 'refAnt'])
    return run(soltab, soltabOut, minFreq, refAnt)


def run( soltab, soltabOut='phasediff', minFreq=0, refAnt='' ):
    """
    Estimate polarization misalignment as delay.

    Parameters
    ----------
    soltabOut : str, optional
        output table name (same solset), by deault "phasediff".

    minFreq : float, optional
        minimum frequency [Hz] to use in estimating the PA. By default, 0 (all freqs).

    refAnt : str, optional
        Reference antenna, by default 'auto'.
    """
    import numpy as np
    from scipy import stats, optimize
    from scipy.ndimage import generic_filter

    logging.info("Finding polarization align for soltab: "+soltab.name)

    solType = soltab.getType()
    if solType != 'phase':
        logging.warning("Soltab type of "+soltab.name+" is of type "+solType+", should be phase. Ignoring.")
        return 1
    
    if refAnt == '': refAnt = 'auto'
    elif refAnt != 'closest' and refAnt != 'auto' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.warning('Reference antenna '+refAnt+' not found. Using: atomatic search.')
        refAnt = 'auto'

    # times and ants needs to be complete or selection is much slower
    times = soltab.getAxisValues('time')

    # create new table
    solset = soltab.getSolset()
    soltabout = solset.makeSoltab(soltype = soltab.getType(), soltabName = soltabOut, axesNames=soltab.getAxesNames(),
                      axesVals=[soltab.getAxisValues(axisName) for axisName in soltab.getAxesNames()],
                      vals=soltab.getValues(retAxesVals = False), weights=soltab.getValues(weight = True, retAxesVals = False))
    soltabout.addHistory('Created by POLALIGN operation from %s.' % soltab.name)

    if 'XX' in soltab.getAxisValues('pol'): pol = 'XX'
    elif 'RR' in soltab.getAxisValues('pol'): pol = 'RR'
    else:
        logging.error('Cannot reference to known polarisation.')
        return 1

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['freq','pol','time'], weight=True, refAnt=refAnt):

        # reorder axes
        vals = reorderAxes( vals, soltab.getAxesNames(), ['pol','freq','time'] )
        weights = reorderAxes( weights, soltab.getAxesNames(), ['pol','freq','time'] )

        if 'RR' in coord['pol'] and 'LL' in coord['pol']:
            coord1 = np.where(coord['pol'] == 'RR')[0][0]
            coord2 = np.where(coord['pol'] == 'LL')[0][0]
        elif 'XX' in coord['pol'] and 'YY' in coord['pol']:
            coord1 = np.where(coord['pol'] == 'XX')[0][0]
            coord2 = np.where(coord['pol'] == 'YY')[0][0]

        if (weights == 0.).all() == True:
            logging.warning('Skipping flagged antenna: '+coord['ant'])
            weights[:] = 0

        elif coord['ant'] == refAnt or (refAnt == 'auto' and coord['ant'] == soltab.cacheAutoRefAnt):
            vals[:] = 0
            weights[:] = 1

        else:
    
            vals[ weights == 0 ] = np.nan # flags (correctly ignored in circmean/circvar or propagated as nan if all flagged, taken care right after usign weights
            if minFreq > 0:
                idx = (coord['freq'] < minFreq)
                vals[:, idx, :] = np.nan

            phase_diff = vals[coord1,:,:] - vals[coord2,:,:]
            phase_diff_mean = stats.circmean(phase_diff, high=1.*np.pi, low=-1.*np.pi, axis=1, nan_policy='omit')
            phase_diff_var = stats.circvar(phase_diff, high=1.*np.pi, low=-1.*np.pi, axis=1, nan_policy='omit')
            phase_diff_var = np.sqrt(phase_diff_var)/np.count_nonzero(~np.isnan(phase_diff), axis=1)**2 # scale the variance by the number of good data points squared
            # remove bad points giving very high weight
            phase_diff_mean[np.isnan(phase_diff_mean)] = 0
            phase_diff_var[ phase_diff_mean == 0 ] = 1e9
            phase_diff_var[ phase_diff_var == 0 ] = 1e9

            # less stable with wraps, but easier code:
            #def func(x, a, b):
            #    return a * x + b
            #popt, pcov = optimize.curve_fit(func, coord['freq'], phase_diff_mean, p0=(1e-9,0.1), sigma=phase_diff_var)

            # more stable with wraps (not really working well):
            #def func(d, freq, y, w):
            #    return np.sum(w*abs(np.cos(d[0]*freq+d[1]) - np.cos(y)) + w*abs(np.sin(d[0]*freq+d[1]) - np.sin(y)))
            #result = optimize.minimize(func, (1e-9, 0.1), args=(coord['freq'], phase_diff_mean, 1/phase_diff_var))
            #popt = result.x

            # more stable with wraps + brute:
            delaycomplex = lambda d, freq, y, w: np.sum(w*abs(np.cos(d[0]*freq+d[1]) - np.cos(y)) + w*abs(np.sin(d[0]*freq+d[1]) - np.sin(y)))
            ranges = (slice(-3e-9, 3e-9, 1e-10),slice(-5,5,0.1))
            popt = optimize.brute(delaycomplex, ranges, finish=optimize.fmin, args=(coord['freq'], phase_diff_mean, 1/phase_diff_var))

            logging.info('%s: average delay: %f ns (offset: %f)' % ( coord['ant'], popt[0]*1e9/(2*np.pi), popt[1]))

            if coord['ant'] == 'UK608HBA':
                print(phase_diff_mean[:100], phase_diff_var[:100], np.count_nonzero(~np.isnan(phase_diff[:100]), axis=1))
                if not 'matplotlib' in sys.modules:
                    import matplotlib as mpl
                    mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
                    mpl.use("Agg")
                    import matplotlib.pyplot as plt

                fig = plt.figure()
                fig.subplots_adjust(wspace=0)
                ax = fig.add_subplot(111)

                # plot rm fit
                plotdelay = lambda delay, offset, freq: np.mod( delay*freq + offset + np.pi, 2.*np.pi) - np.pi
                ax.plot(coord['freq']/1e6, popt[0]*coord['freq'] + popt[1], "-", color='black',  zorder=10, label=r'delay:%f$\nu$ (ns) + %f ' % (popt[0]*1e9/(2*np.pi),popt[1]) )

                ax.plot(coord['freq']/1e6, phase_diff_mean, '.', color='purple' , label='phase difference')

                residual = np.mod(plotdelay(popt[0], popt[1], coord['freq'])-phase_diff_mean + np.pi,2.*np.pi)-np.pi
                ax.plot(coord['freq']/1e6, residual, '.', color='yellow', label='residual')

                ax.set_xlabel('freq [MHz]')
                ax.set_ylabel('phase')
                ax.set_ylim(ymin=-np.pi, ymax=np.pi)

                logging.warning('Save pic: '+coord['ant']+'.png')
                fig.legend()
                plt.savefig(coord['ant']+'.png', bbox_inches='tight')
                del fig

            vals[coord1,:,:] = 0
            phase = np.mod(popt[0]*coord['freq'] + popt[1] + np.pi, 2.*np.pi) - np.pi
            vals[coord2,:,:] = -1.*np.repeat(np.expand_dims(phase, axis=1), repeats=len(coord['time']), axis=1)
            weights[:] = 1

        # reorder axes back to the original order, needed for setValues
        vals = reorderAxes( vals, ['pol','freq','time'], [ax for ax in soltab.getAxesNames() if ax in ['pol','freq','time']] )
        weights = reorderAxes( weights, ['pol','freq','time'], [ax for ax in soltab.getAxesNames() if ax in ['pol','freq','time']] )
        soltabout.setSelection(**coord)
        soltabout.setValues( vals )
        soltabout.setValues( weights, weight=True )

    return 0
