#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Faraday rotation extraction

import logging
from losoto.operations_lib import *
logging.debug('Loading FARADAY module.')

def run( step, parset, H ):
    """
    Separate phase solutions into FR, Clock and TEC.

    The Clock and TEC values are stored in the specified output soltab with type 'clock', 'tec', 'FR'.
    """
    from losoto.h5parm import solFetcher, solWriter
    import numpy as np
    import scipy.optimize

    rmwavcomplex = lambda RM, wav, y: abs(np.cos(2.*RM[0]*wav*wav)  - np.cos(y)) + abs(np.sin(2.*RM[0]*wav*wav)  - np.sin(y))
    c = 2.99792458e8

    # get involved solsets using local step values or global values or all
    soltabs = getParSoltabs( step, parset, H )

    refAnt = parset.getString('.'.join(["LoSoTo.Steps", step, "RefAnt"]), '' )
    max_residual = parset.getFloat('.'.join(["LoSoTo.Steps", step, "MaxResidual"]), 1 )

    for t, soltab in enumerate(openSoltabs( H, soltabs )):
        logging.info("--> Working on soltab: "+soltab._v_name)
        sf = solFetcher(soltab)

        # times and ants needs to be complete or selection is much slower
        times = sf.getAxisValues('time')
        ants = sf.getAxisValues('ant')

        # this will make a selection for the getValues() and getValuesIter()
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        # some checks
        solType = sf.getType()
        if solType != 'phase':
           logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase. Ignoring.")
           continue

        if refAnt != '' and not refAnt in ants:
            logging.error('Reference antenna '+refAnt+' not found. Using: '+ants[1])
            refAnt = ants[1]
        if refAnt == '': refAnt = ants[1]

        if 'XX' in sf.getAxisValues('pol') and 'YY' in sf.getAxisValues('pol'):
            logging.warning('Linear polarization detected, LoSoTo assumes XX->RR and YY->LL.')

        # create new table
        solsetname = soltabs[t].split('/')[0]
        st = H.makeSoltab(solsetname, 'rotationmeasure',
                                 axesNames=['ant','time'], axesVals=[ants, times],
                                 vals=np.zeros((len(ants),len(times))),
                                 weights=np.ones((len(ants),len(times))))
        sw = solWriter(st)
        sw.addHistory('Created by FARADAY operation.')
            
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=['freq','pol','time'], weight=True, reference = refAnt):

            if len(coord['freq']) < 10:
                logging.error('Faraday rotation estimation needs at least 10 frequency channels, preferably distributed over a wide range.')
                return 1

            fitrm = np.zeros(len(times))
            fitweights = np.ones(len(times)) # all unflagged to start
            fitrmguess = 0.001 # good guess

            if 'RR' in coord['pol'] and 'LL' in coord['pol']:
                coord_rr = np.where(coord['pol'] == 'RR')[0][0]
                coord_ll = np.where(coord['pol'] == 'LL')[0][0]
            elif 'XX' in coord['pol'] and 'YY' in coord['pol']:
                coord_rr = np.where(coord['pol'] == 'XX')[0][0]
                coord_ll = np.where(coord['pol'] == 'YY')[0][0]
            else:
                logging.error("Cannot proceed with Faraday estimation with polarizations: "+str(coord['pol']))
                return 1

            if not coord['ant'] == refAnt:
                logging.debug('Working on ant: '+coord['ant']+'...')

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
    
                        if len(freq) < 10:
                            fitweights[t] = 0
                            logging.warning('No valid data found for Faraday fitting for antenna: '+coord['ant']+' at timestamp '+str(t))
                            continue
            
                        if (len(idx) - len(freq))/len(freq) > 1/4.:
                            logging.debug('High number of filtered out data points for the timeslot '+str(t)+': '+str(len(weights[0,:,t]) - len(idx)))

                        # RR-LL to be consistent with BBS/NDPPP
                        phase_diff  = (phase_rr - phase_ll)      # not divide by 2 otherwise jump problem, then later fix this
                        wav = c/freq
        
                        fitresultrm_wav, success = scipy.optimize.leastsq(rmwavcomplex, [fitrmguess], args=(wav, phase_diff))
                        # fractional residual
                        residual = np.mean(np.abs(np.mod((2.*fitresultrm_wav*wav*wav)-phase_diff + np.pi, 2.*np.pi) - np.pi))

    #                    print "t:", t, "result:", fitresultrm_wav, "residual:", residual
    
                        if max_residual == 0 or residual < max_residual:
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
                            print "Plotting"
                            if not 'matplotlib' in sys.modules:
                                import matplotlib as mpl
                                mpl.rc('font',size =8 )
                                mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
                                mpl.use("Agg")
                            import matplotlib.pyplot as plt
    
                            fig = plt.figure()
                            fig.subplots_adjust(wspace=0)
                            ax = fig.add_subplot(110)
    
                            # plot rm fit
                            plotrm = lambda RM, wav: np.mod( (2.*RM*wav*wav) + np.pi, 2.*np.pi) - np.pi # notice the factor of 2
                            ax.plot(freq, plotrm(fitresultrm_wav, c/freq[:]), "-", color='purple')
    
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

            sw.setSelection(ant=coord['ant'], time=coord['time'])
            sw.setValues( np.expand_dims(fitrm, axis=1) )
            sw.setValues( np.expand_dims(fitweights, axis=1), weight=True )

        del st
        del sw        
        del sf
    return 0
