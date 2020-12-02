#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading FARADAY module.')

def _run_parser(soltab, parser, step):
    soltabOut = parser.getstr( step, 'soltabOut', 'rotationmeasure000' )
    refAnt = parser.getstr( step, 'refAnt', '')
    maxResidual = parser.getfloat( step, 'maxResidual', 1. )

    parser.checkSpelling( step, soltab, ['soltabOut', 'refAnt', 'maxResidual'])
    return run(soltab, soltabOut, refAnt, maxResidual)


def run( soltab, soltabOut='rotationmeasure000', refAnt='', maxResidual=1. ):
    """
    Faraday rotation extraction from either a rotation table or a circular phase (of which the operation get the polarisation difference).

    Parameters
    ----------
    
    soltabOut : str, optional
        output table name (same solset), by deault "rotationmeasure000".
        
    refAnt : str, optional
        Reference antenna, by default the first.

    maxResidual : float, optional
        Max average residual in radians before flagging datapoint, by default 1. If 0: no check.

    """
    import numpy as np
    import scipy.optimize

    rmwavcomplex = lambda RM, wav, y: abs(np.cos(2.*RM[0]*wav*wav)  - np.cos(y)) + abs(np.sin(2.*RM[0]*wav*wav)  - np.sin(y))
    c = 2.99792458e8

    logging.info("Find FR for soltab: "+soltab.name)

    # input check
    solType = soltab.getType()
    if solType == 'phase':
        returnAxes = ['pol','freq','time']
        if 'RR' in soltab.getAxisValues('pol') and 'LL' in soltab.getAxisValues('pol'):
            coord_rr = np.where(soltab.getAxisValues('pol') == 'RR')[0][0]
            coord_ll = np.where(soltab.getAxisValues('pol') == 'LL')[0][0]
        elif 'XX' in soltab.getAxisValues('pol') and 'YY' in soltab.getAxisValues('pol'):
            logging.warning('Linear polarization detected, LoSoTo assumes XX->RR and YY->LL.')
            coord_rr = np.where(soltab.getAxisValues('pol') == 'XX')[0][0]
            coord_ll = np.where(soltab.getAxisValues('pol') == 'YY')[0][0]
        else:
            logging.error("Cannot proceed with Faraday estimation with polarizations: "+str(coord['pol']))
            return 1
    elif solType == 'rotation':
        returnAxes = ['freq','time']
    else:
       logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase or rotation. Ignoring.")
       return 1

    ants = soltab.getAxisValues('ant')
    if refAnt != '' and refAnt != 'closest' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.warning('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[0]
    if refAnt == '': refAnt = ants[0]

    # times and ants needs to be complete or selection is much slower
    times = soltab.getAxisValues('time')

    # create new table
    solset = soltab.getSolset()
    soltabout = solset.makeSoltab('rotationmeasure', soltabName = soltabOut,
                             axesNames=['ant','time'], axesVals=[ants, times],
                             vals=np.zeros((len(ants),len(times))),
                             weights=np.ones((len(ants),len(times))))
    soltabout.addHistory('Created by FARADAY operation from %s.' % soltab.name)

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=returnAxes, weight=True, refAnt=refAnt):

        if len(coord['freq']) < 10:
            logging.error('Faraday rotation estimation needs at least 10 frequency channels, preferably distributed over a wide range.')
            return 1

        # reorder axes
        vals = reorderAxes( vals, soltab.getAxesNames(), returnAxes )
        weights = reorderAxes( weights, soltab.getAxesNames(), returnAxes )
        weights[np.isnan(vals)] = 0.

        fitrm = np.zeros(len(times))
        fitweights = np.ones(len(times)) # all unflagged to start
        fitrmguess = 0.001 # good guess

        if not coord['ant'] == refAnt:
            logging.debug('Working on ant: '+coord['ant']+'...')

            if (weights == 0.).all() == True:
                logging.warning('Skipping flagged antenna: '+coord['ant'])
                fitweights[:] = 0
            else:

                for t, time in enumerate(times):

                    if solType == 'phase':
                        idx       = ((weights[coord_rr,:,t] != 0.) & (weights[coord_ll,:,t] != 0.))
                        freq      = np.copy(coord['freq'])[idx]
                        phase_rr  = vals[coord_rr,:,t][idx]
                        phase_ll  = vals[coord_ll,:,t][idx]
                        # RR-LL to be consistent with BBS/NDPPP
                        phase_diff  = (phase_rr - phase_ll)      # not divide by 2 otherwise jump problem, then later fix this
                    else: # rotation table
                        idx        = ((weights[:,t] != 0.) & (weights[:,t] != 0.))
                        freq       = np.copy(coord['freq'])[idx]
                        phase_diff = 2.*vals[:,t][idx] # a rotation is between -pi and +pi

                    if len(freq) < 20:
                        fitweights[t] = 0
                        logging.warning('No valid data found for Faraday fitting for antenna: '+coord['ant']+' at timestamp '+str(t))
                        continue
        
                    # if more than 1/4 of chans are flagged
                    if (len(idx) - len(freq))/float(len(idx)) > 1/4.:
                        logging.debug('High number of filtered out data points for the timeslot %i: %i/%i' % (t, len(idx) - len(freq), len(idx)) )

                    wav = c/freq
    
                    fitresultrm_wav, success = scipy.optimize.leastsq(rmwavcomplex, [fitrmguess], args=(wav, phase_diff))
                    # fractional residual
                    residual = np.nanmean(np.abs(np.mod((2.*fitresultrm_wav*wav*wav)-phase_diff + np.pi, 2.*np.pi) - np.pi))

#                    print "t:", t, "result:", fitresultrm_wav, "residual:", residual

                    if maxResidual == 0 or residual < maxResidual:
                        fitrmguess = fitresultrm_wav[0]
                        weight = 1
                    else:       
                        # high residual, flag
                        logging.warning('Bad solution for ant: '+coord['ant']+' (time: '+str(t)+', resdiaul: '+str(residual)+').')
                        weight = 0

                    fitrm[t] = fitresultrm_wav[0]
                    fitweights[t] = weight

                    # Debug plot
                    doplot = False
                    if doplot and coord['ant'] == 'RS310LBA' and t%10==0:
                        print("Plotting")
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
                        plotrm = lambda RM, wav: np.mod( (2.*RM*wav*wav) + np.pi, 2.*np.pi) - np.pi # notice the factor of 2
                        ax.plot(freq, plotrm(fitresultrm_wav, c/freq[:]), "-", color='purple')

                        if solType == 'phase':
                            ax.plot(freq, np.mod(phase_rr + np.pi, 2.*np.pi) - np.pi, 'ob' )
                            ax.plot(freq, np.mod(phase_ll + np.pi, 2.*np.pi) - np.pi, 'og' )
                        ax.plot(freq, np.mod(phase_diff + np.pi, 2.*np.pi) - np.pi , '.', color='purple' )                           
     
                        residual = np.mod(plotrm(fitresultrm_wav, c/freq[:])-phase_diff+np.pi,2.*np.pi)-np.pi
                        ax.plot(freq, residual, '.', color='yellow')
        
                        ax.set_xlabel('freq')
                        ax.set_ylabel('phase')
                        ax.set_ylim(ymin=-np.pi, ymax=np.pi)
    
                        logging.warning('Save pic: '+str(t)+'_'+coord['ant']+'.png')
                        plt.savefig(str(t)+'_'+coord['ant']+'.png', bbox_inches='tight')
                        del fig

        soltabout.setSelection(ant=coord['ant'], time=coord['time'])
        soltabout.setValues( np.expand_dims(fitrm, axis=1) )
        soltabout.setValues( np.expand_dims(fitweights, axis=1), weight=True )

    return 0
