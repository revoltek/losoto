#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.operations_lib import *

logging.debug('Loading DELAY module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'delay' )
    refAnt = parser.getstr( step, 'refAnt', '')
    maxResidual = parser.getfloat( step, 'maxResidual', 1. )
    return run(soltab, soltabOut, refAnt, maxResidual)


def run( soltab, soltabOut='delay', refAnt='', maxResidual=1. ):
    """
    Differential delay extraction.

    Parameters
    ----------
    soltabOut : str, optional
        output table name (same solset), by deault "delay".

    refAnt : str, optional
        Reference antenna, by default the first.

    maxResidual : float, optional
        Max average residual in radians before flagging datapoint, by default 1. If 0: no check.

    """
    import numpy as np
    import scipy.optimize

    dcomplex = lambda d, freq, y: abs(np.cos(d[0]*freq)  - np.cos(y)) + abs(np.sin(d[0]*freq)  - np.sin(y))
    #dcomplex2 = lambda d, freq, y: abs(np.cos(d[0]*freq + d[1])  - np.cos(y)) + abs(np.sin(d[0]*freq + d[1])  - np.sin(y))

    logging.info("Find DELAY for soltab: "+soltab.name)

    # input check
    solType = soltab.getType()
    if solType != 'phase':
       logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase. Ignoring.")
       return 1

    ants = soltab.getAxisValues('ant')
    if refAnt != '' and not refAnt in ants:
        logging.error('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[0]
    if refAnt == '': refAnt = ants[0]

    # times and ants needs to be complete or selection is much slower
    times = soltab.getAxisValues('time')

    # create new table
    solset = soltab.getSolset()
    soltabout = solset.makeSoltab(soltype = soltab.getType(), soltabName = soltabOut, axesNames=soltab.getAxesNames(), \
                      axesVals=[soltab.getAxisValues(axisName) for axisName in soltab.getAxesNames()], \
                      vals=soltab.getValues(retAxesVals = False), weights=soltab.getValues(weight = True, retAxesVals = False))
    soltabout.addHistory('Created by DELAY operation.')
        
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['freq','pol','time'], weight=True, reference=refAnt):

        if len(coord['freq']) < 10:
            logging.error('Delay estimation needs at least 10 frequency channels, preferably distributed over a wide range.')
            return 1

        # reorder axes
        vals = reorderAxes( vals, soltab.getAxesNames(), ['pol','freq','time'] )
        weights = reorderAxes( weights, soltab.getAxesNames(), ['pol','freq','time'] )

        fitd = np.zeros(len(times))
        fitweights = np.ones(len(times)) # all unflagged to start
        fitdguess = 1e-9 # good guess

        if 'RR' in coord['pol'] and 'LL' in coord['pol']:
            coord_rr = np.where(coord['pol'] == 'RR')[0][0]
            coord_ll = np.where(coord['pol'] == 'LL')[0][0]
        elif 'XX' in coord['pol'] and 'YY' in coord['pol']:
            coord_rr = np.where(coord['pol'] == 'XX')[0][0]
            coord_ll = np.where(coord['pol'] == 'YY')[0][0]
        else:
            logging.error("Cannot proceed with delay estimation with polarizations: "+str(coord['pol']))
            return 1

        if not coord['ant'] == refAnt:

            if (weights == 0.).all() == True:
                logging.warning('Skipping flagged antenna: '+coord['ant'])
                fitweights[:] = 0
            else:

                for t, time in enumerate(times):

                    # apply flags
                    idx       = ((weights[coord_rr,:,t] != 0.) & (weights[coord_ll,:,t] != 0.))
                    freq      = np.copy(coord['freq'])[idx]
                    phase_rr  = vals[coord_rr,:,t][idx]
                    phase_ll  = vals[coord_ll,:,t][idx]

                    if len(freq) < 30:
                        fitweights[t] = 0
                        logging.warning('No valid data found for delay fitting for antenna: '+coord['ant']+' at timestamp '+str(t))
                        continue
        
                    # if more than 1/4 of chans are flagged
                    if (len(idx) - len(freq))/float(len(idx)) > 1/4.:
                        logging.debug('High number of filtered out data points for the timeslot %i: %i/%i' % (t, len(idx) - len(freq), len(idx)) )

                    # RR-LL to be consistent with BBS/NDPPP
                    phase_diff  = (phase_rr - phase_ll)      # not divide by 2 otherwise jump problem, then later fix this
                    phase_diff = np.unwrap(phase_diff)
    
                    #fitresultd2, success = scipy.optimize.leastsq(dcomplex2, [fitdguess,-1.], args=(freq, phase_diff))
                    #numjumps = np.around(fitresultd2[1]/(2*np.pi))
                    #if numjumps > 0: print fitresultd, numjumps
                    #phase_diff += numjumps * 2*np.pi

                    best_residual = np.inf
                    for jump in [-2,-1,0,1,2]:
                        fitresultd, success = scipy.optimize.leastsq(dcomplex, [fitdguess], args=(freq, phase_diff - jump * 2*np.pi))
                        # fractional residual
                        residual = np.nanmean(np.abs( (fitresultd[0]*freq) - phase_diff - jump * 2*np.pi ) )
                        if residual < best_residual:
                            best_residual = residual
                            fitd[t] = fitresultd[0]
                            best_jump = jump

                    if maxResidual == 0 or best_residual < maxResidual:
                        fitweights[t] = 1
                        fitdguess = fitresultd[0]
                    else:       
                        # high residual, flag
                        logging.warning('Bad solution for ant: '+coord['ant']+' (time: '+str(t)+', resdiaul: '+str(best_residual)+').')
                        fitweights[t] = 0

                    # Debug plot
                    doplot = True
                    if doplot and (coord['ant'] == 'RS310LBA' or coord['ant'] == 'RS210LBA') and t%100==0:
                        print "Plotting"
                        if not 'matplotlib' in sys.modules:
                            import matplotlib as mpl
                            mpl.rc('font',size =8 )
                            mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
                            mpl.use("Agg")
                        import matplotlib.pyplot as plt

                        fig = plt.figure()
                        fig.subplots_adjust(wspace=0)
                        ax = fig.add_subplot(111)

                        # plot rm fit
                        plotd = lambda d, freq: d*freq # notice the factor of 2
                        ax.plot(freq, plotd(fitd[t], freq[:]), "-", color='purple')

                        ax.plot(freq, np.mod(phase_rr + np.pi, 2.*np.pi) - np.pi, 'ob' )
                        ax.plot(freq, np.mod(phase_ll + np.pi, 2.*np.pi) - np.pi, 'og' )
                        ax.plot(freq, phase_diff - best_jump * 2*np.pi , '.', color='purple' )                           
     
                        residual = np.mod(plotd(fitd[t], freq[:])-phase_diff+np.pi,2.*np.pi)-np.pi
                        ax.plot(freq, residual, '.', color='yellow')
        
                        ax.set_xlabel('freq')
                        ax.set_ylabel('phase')
                        #ax.set_ylim(ymin=-np.pi, ymax=np.pi)
    
                        logging.warning('Save pic: '+str(t)+'_'+coord['ant']+'.png')
                        plt.savefig(str(t)+'_'+coord['ant']+'.png', bbox_inches='tight')
                        del fig

                logging.debug('%s: average delay: %f ns' % (coord['ant'], np.mean(2*fitd)*1e9))
                for t, time in enumerate(times):
                    #vals[:,:,t] = 0.
                    #vals[coord1,:,t] = fitd[t]*np.array(coord['freq'])/2.
                    vals[coord_rr,:,t] = 0
                    vals[coord_ll,:,t] = -1.*(fitd[t]*np.array(coord['freq']))#/2.
                    #weights[:,:,t] = 0.
                    weights[coord_rr,:,t] = fitweights[t]
                    weights[coord_ll,:,t] = fitweights[t]

        # reorder axes back to the original order, needed for setValues
        vals = reorderAxes( vals, ['pol','freq','time'], [ax for ax in soltab.getAxesNames() if ax in ['pol','freq','time']] )
        weights = reorderAxes( weights, ['pol','freq','time'], [ax for ax in soltab.getAxesNames() if ax in ['pol','freq','time']] )
        soltabout.setSelection(**coord)
        soltabout.setValues( vals )
        soltabout.setValues( weights, weight=True )

    return 0
