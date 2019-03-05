#!/usr/bin/env python
# -*- coding: utf-8 -*-


import logging
from losoto.lib_operations import *

logging.debug('Loading TEC module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'tec000' )
    refAnt = parser.getstr( step, 'refAnt', '')
    maxResidual = parser.getfloat( step, 'maxResidual', 2.5 )

    parser.checkSpelling( step, soltab, ['soltabOut', 'refAnt', 'maxResidual'])
    return run(soltab, soltabOut, refAnt, maxResidual)


def run( soltab, soltabOut='tec000', refAnt='', maxResidual=1. ):
    """
    Bruteforce TEC extraction from phase solutions.

    Parameters
    ----------
    soltabOut : str, optional
        output table name (same solset), by deault "tec".

    refAnt : str, optional
        Reference antenna, by default the first.

    maxResidual : float, optional
        Max average residual in radians before flagging datapoint, by default 1. If 0: no check.

    """
    import numpy as np
    import scipy.optimize
    from losoto.lib_unwrap import unwrap_2d

    def mod(d):
        return np.mod(d + np.pi, 2.*np.pi) - np.pi

    drealbrute = lambda d, freq, y: np.sum(np.abs(mod(-8.44797245e9*d/freq) - y)) # bruteforce
    dreal = lambda d, freq, y: mod(-8.44797245e9*d[0]/freq) - y
    #dreal2 = lambda d, freq, y: mod(-8.44797245e9*d[0]/freq + d[1]) - y
    #dcomplex = lambda d, freq, y:  np.sum( ( np.cos(-8.44797245e9*d/freq)  - np.cos(y) )**2 ) +  np.sum( ( np.sin(-8.44797245e9*d/freq)  - np.sin(y) )**2 ) 
    #def dcomplex( d, freq, y, y_pre, y_post):  
    #    return np.sum( ( np.absolute( np.exp(-1j*8.44797245e9*d/freq)  - np.exp(1j*y) ) )**2 ) + \
    #           .5*np.sum( ( np.absolute( np.exp(-1j*8.44797245e9*d/freq)  - np.exp(1j*y_pre) ) )**2 ) + \
    #           .5*np.sum( ( np.absolute( np.exp(-1j*8.44797245e9*d/freq)  - np.exp(1j*y_post) ) )**2 )
    #dcomplex2 = lambda d, freq, y:  abs(np.cos(-8.44797245e9*d[0]/freq + d[1])  - np.cos(y)) + abs(np.sin(-8.44797245e9*d[0]/freq + d[1])  - np.sin(y))

    logging.info("Find TEC for soltab: "+soltab.name)

    # input check
    solType = soltab.getType()
    if solType != 'phase':
       logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase. Ignoring.")
       return 1

    ants = soltab.getAxisValues('ant')
    if refAnt != '' and refAnt != 'closest' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.error('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[0]
    if refAnt == '': refAnt = ants[0]

    # times and ants needs to be complete or selection is much slower
    times = soltab.getAxisValues('time')

    # create new table
    solset = soltab.getSolset()
    soltabout = solset.makeSoltab(soltype = 'tec', soltabName = soltabOut, axesNames=['ant','time'], \
                      axesVals=[soltab.getAxisValues(axisName) for axisName in ['ant','time']], \
                      vals=np.zeros(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('time'))), \
                      weights=np.ones(shape=(soltab.getAxisLen('ant'),soltab.getAxisLen('time'))) )
    soltabout.addHistory('Created by TEC operation from %s.' % soltab.name)
        
    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['freq','time'], weight=True, reference=refAnt):

        if len(coord['freq']) < 10:
            logging.error('Delay estimation needs at least 10 frequency channels, preferably distributed over a wide range.')
            return 1

        # reorder axes
        vals = reorderAxes( vals, soltab.getAxesNames(), ['freq','time'] )
        weights = reorderAxes( weights, soltab.getAxesNames(), ['freq','time'] )

        fitd = np.zeros(len(times))
        fitweights = np.ones(len(times)) # all unflagged to start
        #fitdguess = 0.01 # good guess
        ranges = (-0.5,0.5)
        Ns = 1000

        if not coord['ant'] == refAnt:

            if (weights == 0.).all() == True:
                logging.warning('Skipping flagged antenna: '+coord['ant'])
                fitweights[:] = 0
            else:

                # unwrap 2d timexfreq
                #flags = np.array((weights == 0), dtype=bool)
                #vals = unwrap_2d(vals, flags, coord['freq'], coord['time'])

                for t, time in enumerate(times):

                    # apply flags
                    idx       = (weights[:,t] != 0.)
                    freq      = np.copy(coord['freq'])[idx]

                    if t == 0: phaseComb_pre  = vals[idx,0]
                    else: phaseComb_pre  = vals[idx,t-1]
                    if t == len(times)-1: phaseComb_post  = vals[idx,-1]
                    else: phaseComb_post  = vals[idx,t+1]
                    
                    phaseComb  = vals[idx,t]

                    if len(freq) < 10:
                        fitweights[t] = 0
                        logging.warning('No valid data found for delay fitting for antenna: '+coord['ant']+' at timestamp '+str(t))
                        continue
        
                    # if more than 1/4 of chans are flagged
                    if (len(idx) - len(freq))/float(len(idx)) > 1/4.:
                        logging.debug('High number of filtered out data points for the timeslot %i: %i/%i' % (t, len(idx) - len(freq), len(idx)) )

                    # least square 2
                    #fitresultd2, success = scipy.optimize.leastsq(dreal2, [fitdguess,0.], args=(freq, phaseComb))
                    #numjumps = np.around(fitresultd2[1]/(2*np.pi))
                    #print 'best jumps:', numjumps
                    #phaseComb -= numjumps * 2*np.pi
                    #fitresultd, success = scipy.optimize.leastsq(dreal, [fitresultd2[0]], args=(freq, phaseComb))

                    # least square 1
                    #fitresultd, success = scipy.optimize.leastsq(dreal, [fitdguess], args=(freq, phaseComb))

                    # hopper
                    #fitresultd = scipy.optimize.basinhopping(dreal, [fitdguess], T=1, minimizer_kwargs={'args':(freq, phaseComb)})
                    #fitresultd = [fitresultd.x]

                    #best_residual = np.nanmean(np.abs( mod(-8.44797245e9*fitresultd[0]/freq) - phaseComb ) )

                    #best_residual = np.inf
                    #for jump in [-2,-1,0,1,2]:
                    #    fitresultd, success = scipy.optimize.leastsq(dreal, [fitdguess], args=(freq, phaseComb - jump * 2*np.pi))
                    #    print fitresultd
                    #    # fractional residual
                    #    residual = np.nanmean(np.abs( (-8.44797245e9*fitresultd[0]/freq) - phaseComb - jump * 2*np.pi ) )
                    #    if residual < best_residual:
                    #        best_residual = residual
                    #        fitd[t] = fitresultd[0]
                    #        best_jump = jump

                    # brute force
                    fitresultd = scipy.optimize.brute(drealbrute, ranges=(ranges,), Ns=Ns, args=(freq, phaseComb))
                    fitresultd, success = scipy.optimize.leastsq(dreal, fitresultd, args=(freq, phaseComb))
                    best_residual = np.nanmean(np.abs( mod(-8.44797245e9*fitresultd[0]/freq) - phaseComb ) )

                    fitd[t] = fitresultd[0]
                    if maxResidual == 0 or best_residual < maxResidual:
                        fitweights[t] = 1
                        #fitdguess = fitresultd[0]
                        ranges = (fitresultd[0]-0.05,fitresultd[0]+0.05)
                        Ns = 100
                    else:       
                        # high residual, flag and reset initial guess
                        logging.warning('Bad solution for ant: '+coord['ant']+' (time: '+str(t)+', resdiaul: '+str(best_residual)+').')
                        fitweights[t] = 0
                        #fitdguess = 0.01
                        ranges = (-0.5,0.5)
                        Ns = 1000


                    # Debug plot
                    doplot = False
                    if doplot and (coord['ant'] == 'RS509LBA' or coord['ant'] == 'RS210LBA') and t%50==0:
                        print("Plotting")
                        if not 'matplotlib' in sys.modules:
                            import matplotlib as mpl
                            mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
                            mpl.use("Agg")
                        import matplotlib.pyplot as plt

                        fig = plt.figure()
                        fig.subplots_adjust(wspace=0)
                        ax = fig.add_subplot(111)

                        # plot rm fit
                        plotd = lambda d, freq: -8.44797245e9*d/freq
                        ax.plot(freq, plotd(fitresultd[0], freq[:]), "-", color='purple')
                        ax.plot(freq, mod(plotd(fitresultd[0], freq[:])), ":", color='purple')

                        #ax.plot(freq, vals[idx,t], '.b' )
                        #ax.plot(freq, phaseComb + numjumps * 2*np.pi, 'x', color='purple' )                           
                        ax.plot(freq, phaseComb, 'o', color='purple' )                           
     
                        residual = mod( plotd(fitd[t], freq[:]) - phaseComb)
                        ax.plot(freq, residual, '.', color='orange')
        
                        ax.set_xlabel('freq')
                        ax.set_ylabel('phase')
                        #ax.set_ylim(ymin=-np.pi, ymax=np.pi)
    
                        logging.warning('Save pic: '+str(t)+'_'+coord['ant']+'.png')
                        plt.savefig(str(t)+'_'+coord['ant']+'.png', bbox_inches='tight')
                        del fig

                logging.info('%s: average tec: %f TECU' % (coord['ant'], np.mean(2*fitd)))

        # reorder axes back to the original order, needed for setValues
        soltabout.setSelection(ant=coord['ant'])
        soltabout.setValues( fitd )
        soltabout.setValues( fitweights, weight=True )

    return 0
