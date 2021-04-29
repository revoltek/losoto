#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading POLALIGN module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'phasediff' )
    maxResidual = parser.getfloat( step, 'maxResidual', 1. )
    fitOffset = parser.getbool( step, 'fitOffset', False )
    average = parser.getbool( step, 'average', False )
    replace = parser.getbool( step, 'replace', False )
    minFreq = parser.getfloat( step, 'minFreq', 0 )
    refAnt = parser.getstr( step, 'refAnt', '' )

    parser.checkSpelling( step, soltab, ['soltabOut', 'maxResidual', 'fitOffset', 'average', 'replace', 'minFreq', 'refAnt'])
    return run(soltab, soltabOut, maxResidual, fitOffset, average, replace, minFreq, refAnt)


def run( soltab, soltabOut='phasediff', maxResidual=1., fitOffset=False, average=False, replace=False, minFreq=0, refAnt='' ):
    """
    Estimate polarization misalignment as delay.

    Parameters
    ----------
    soltabOut : str, optional
        output table name (same solset), by deault "phasediff".

    maxResidual : float, optional
        Maximum acceptable rms of the residuals in radians before flagging, by default 1. If 0: No check.

    fitOffset : bool, optional
        Assume that together with a delay each station has also a differential phase offset (important for old LBA observations). By default False.    

    average : bool, optional
        Mean-average in time the resulting delays/offset. By default False.    
    
    replace : bool, optional
        replace using smoothed value instead of flag bad data? Smooth must be active. By default, False.

    minFreq : float, optional
        minimum frequency [Hz] to use in estimating the PA. By default, 0 (all freqs).

    refAnt : str, optional
        Reference antenna, by default the first.
    """
    import numpy as np
    import scipy.optimize
    from scipy import stats
    from scipy.ndimage import generic_filter

    logging.info("Finding polarization align for soltab: "+soltab.name)

    delaycomplex = lambda d, freq, y: abs(np.cos(d[0]*freq)  - np.cos(y)) + abs(np.sin(d[0]*freq)  - np.sin(y))
    #delaycomplex = lambda d, freq, y: abs(d[0]*freq  - y)

    solType = soltab.getType()
    if solType != 'phase':
        logging.warning("Soltab type of "+soltab.name+" is of type "+solType+", should be phase. Ignoring.")
        return 1
    
    if refAnt != '' and refAnt != 'closest' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.warning('Reference antenna '+refAnt+' not found. Using: '+soltab.getAxisValues('ant')[1])
        refAnt = soltab.getAxisValues('ant')[1]
    if refAnt == '': refAnt = soltab.getAxisValues('ant')[1]

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
        else:

            fit_delays=[]; fit_offset=[]; fit_weights=[]
            for t, time in enumerate(times):

                # apply flags
                idx       = ( (weights[coord1,:,t] != 0.) & (weights[coord2,:,t] != 0.) & (coord['freq'] > minFreq) )
                freq      = np.copy(coord['freq'])[idx]
                phase1    = vals[coord1,:,t][idx]
                phase2    = vals[coord2,:,t][idx]

                if len(freq) < 30:
                    fit_weights.append(0.)
                    fit_delays.append(0.)
                    fit_offset.append(0.)
                    logging.debug('Not enough unflagged point for the timeslot '+str(t))
                    continue
    
                # if more than 1/2 of chans are flagged
                if (len(idx) - len(freq))/float(len(idx)) > 1/2.:
                    logging.debug('High number of filtered out data points for the timeslot %i: %i/%i' % (t, len(idx) - len(freq), len(idx)) )
    
                phase_diff = phase1 - phase2
                phase_diff = np.mod(phase_diff + np.pi, 2.*np.pi) - np.pi
                phase_diff = np.unwrap(phase_diff)

                A = np.vstack([freq, np.ones(len(freq))]).T
                fitresultdelay = np.linalg.lstsq(A, phase_diff.T)[0]
                # get the closest n*(2pi) to the intercept and refit with only 1 parameter
                if not fitOffset:
                    numjumps = np.around(fitresultdelay[1]/(2*np.pi))
                    A = np.reshape(freq, (-1,1)) # no b
                    phase_diff -= numjumps * 2 * np.pi
                    fitresultdelay = np.linalg.lstsq(A, phase_diff.T)[0]
                    fitresultdelay = [fitresultdelay[0],0.] # set offset to 0 to keep the rest of the script equal

                # fractional residual
                residual = np.mean(np.abs( fitresultdelay[0]*freq + fitresultdelay[1] - phase_diff ))

                fit_delays.append(fitresultdelay[0])
                fit_offset.append(fitresultdelay[1])
                if maxResidual == 0 or residual < maxResidual:
                    fit_weights.append(1.)
                else:       
                    # high residual, flag
                    logging.debug('Bad solution for ant: '+coord['ant']+' (time: '+str(t)+', residual: '+str(residual)+') -> ignoring.')
                    fit_weights.append(0.)

                # Debug plot
                doplot = False
                if doplot and t%100==0 and (coord['ant'] == 'RS310LBA' or coord['ant'] == 'CS301LBA'):
                # if doplot and t%10==0 and (coord['ant'] == 'W04'):
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
                    ax.plot(freq, fitresultdelay[0]*freq + fitresultdelay[1], "-", color='black',  zorder=10, label=r'delay:%f$\nu$ (ns) + %f ' % (fitresultdelay[0]*1e9,fitresultdelay[1]) )

                    ax.plot(freq, np.mod(phase1 + np.pi, 2.*np.pi) - np.pi, 'ob' , label='phase XX/RR')
                    ax.plot(freq, np.mod(phase2 + np.pi, 2.*np.pi) - np.pi, 'og' , label='phase YY/LL' )
                    #ax.plot(freq, np.mod(phase_diff + np.pi, 2.*np.pi) - np.pi, '.', color='purple' )                           
                    ax.plot(freq, phase_diff, '.', color='purple' , label='phase difference')
 
                    residual = np.mod(plotdelay(fitresultdelay[0], fitresultdelay[1], freq)-phase_diff + np.pi,2.*np.pi)-np.pi
                    ax.plot(freq, residual, '.', color='yellow', label='residual')
    
                    ax.set_xlabel('freq')
                    ax.set_ylabel('phase')
                    #ax.set_ylim(ymin=-np.pi, ymax=np.pi)

                    logging.warning('Save pic: '+str(t)+'_'+coord['ant']+'.png')
                    fig.legend()
                    plt.savefig(coord['ant']+'_'+str(t)+'.png', bbox_inches='tight')
                    del fig
            # end cycle in time

            fit_weights = np.array(fit_weights)
            fit_delays = np.array(fit_delays)
            fit_offset = np.array(fit_offset)

            # avg in time
            if average:
                fit_delays_bkp = fit_delays[ fit_weights == 0 ]
                fit_offset_bkp = fit_offset[ fit_weights == 0 ]
                np.putmask(fit_delays, fit_weights == 0, np.nan)
                np.putmask(fit_offset, fit_weights == 0, np.nan)
                fit_delays[:] = np.nanmean(fit_delays)
                # angle mean
                fit_offset[:] = np.angle( np.nansum( np.exp(1j*fit_offset) ) / np.count_nonzero(~np.isnan(fit_offset)) )

                if replace:
                    fit_weights[ fit_weights == 0 ] = 1.
                    fit_weights[ np.isnan(fit_delays) ] = 0. # all the size was flagged cannot estrapolate value
                else:
                    fit_delays[ fit_weights == 0 ] = fit_delays_bkp
                    fit_offset[ fit_weights == 0 ] = fit_offset_bkp

            logging.info('%s: average delay: %f ns (offset: %f)' % ( coord['ant'], np.mean(fit_delays)*1e9, np.mean(fit_offset)))
            for t, time in enumerate(times):
                #vals[:,:,t] = 0.
                #vals[coord1,:,t] = fit_delays[t]*np.array(coord['freq'])/2.
                vals[coord1,:,t] = 0
                phase = np.mod(fit_delays[t]*coord['freq'] + fit_offset[t] + np.pi, 2.*np.pi) - np.pi
                vals[coord2,:,t] = -1.*phase#/2.
                #weights[:,:,t] = 0.
                weights[coord1,:,t] = fit_weights[t]
                weights[coord2,:,t] = fit_weights[t]

        # reorder axes back to the original order, needed for setValues
        vals = reorderAxes( vals, ['pol','freq','time'], [ax for ax in soltab.getAxesNames() if ax in ['pol','freq','time']] )
        weights = reorderAxes( weights, ['pol','freq','time'], [ax for ax in soltab.getAxesNames() if ax in ['pol','freq','time']] )
        soltabout.setSelection(**coord)
        soltabout.setValues( vals )
        soltabout.setValues( weights, weight=True )

    return 0
