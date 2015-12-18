#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Clock/tec/FR separation module

import logging
from losoto.operations_lib import *
logging.debug('Loading CLOCKTEC2 module.')

c = 2.99792458e8

#def fit_dFR(phases_rr, phases_ll, weights_rr, weights_ll, freq):

def initGuess(sol, ants, freqs, phase, refine=1):
        import numpy as np

        # FIND INTIAL GUESS
        if sol == None:
            sol = np.zeros((len(ants), 2))
            for a, ant in enumerate(ants):
                chis = []; dTECs = []; dclocks = []
                if 'CS' in ant:
                    dclockrange = np.arange(-10e-9, 10e-9, 1.e-9/refine)
                    dTECrange = np.arange(-0.01, 0.01, 0.001/refine)
                else: 
                    dclockrange = np.arange(-200e-9, 200e-9, 5.e-9/refine)
                    dTECrange = np.arange(-0.5, 0.5, 0.01/refine)

                for dTEC in dTECrange:
                    for dclock in dclockrange:
                        phase_model = np.mod( (4.*np.pi*dclock*freqs) - (2.*8.44797245e9*dTEC/freqs), 2*np.pi )  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
                        phase_data  = np.mod( phase[a,:], 2*np.pi )
                        angle       = np.pi - np.abs(np.abs(phase_model - phase_data) - np.pi)
                        chis.append(np.sum(angle))
                        dTECs.append(dTEC)
                        dclocks.append(dclock)
    
                idx = chis.index(min(chis))
                sol[a,:] = [dTECs[idx], dclocks[idx]]
    
        for a, ant in enumerate(ants):
            chis = []; dTECs = []; dclocks = []
            if 'CS' in ant:
                dclockrange = np.arange(sol[a,1]-1.5e-9, sol[a,1]+1.5e-9, 0.1e-9/refine)
                dTECrange = np.arange(sol[a,0]-0.0015, sol[a,0]+0.0015, 0.0001/refine)
            else: 
                dclockrange = np.arange(sol[a,1]-7.5e-9, sol[a,1]+7.5e-9, 1.e-9/refine)
                dTECrange = np.arange(sol[a,0]-0.015, sol[a,0]+0.015, 0.001/refine)

            for dTEC in dTECrange:
                for dclock in dclockrange:
                    phase_model = np.mod( (4.*np.pi*dclock*freqs) - (2.*8.44797245e9*dTEC/freqs), 2*np.pi )  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
                    phase_data  = np.mod( phase[a,:], 2*np.pi )
                    angle       = np.pi - np.abs(np.abs(phase_model - phase_data) - np.pi)
                    chis.append(np.sum(angle))
                    dTECs.append(dTEC)
                    dclocks.append(dclock)
    
            idx = chis.index(min(chis))
            sol[a,:] = [dTECs[idx], dclocks[idx]]

        return sol
 

def fit_dTEC_dclock_dFR(phases, weights, coord):
        import numpy as np
        from lofar.expion import baselinefitting  as fitting
    
        # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
        par3complex = lambda p, freq, y: abs(np.cos((4.*np.pi*p[1]*freq) - (2.*8.44797245e9*p[0]/freq) + p[2]) - np.cos(y)) + abs(np.sin((4.*np.pi*p[1]*freq) - (2.*8.44797245e9*p[0]/freq) + p[2]) - np.sin(y))
        rmwavcomplex = lambda RM, wav, y: abs(np.cos(2.*RM[0]*wav*wav)  - np.cos(y)) + abs(np.sin(2.*RM[0]*wav*wav)  - np.sin(y))

        times = coord['time']
        ants = coord['ant']
        freqs = coord['freq']
        fittec = np.zeros((len(ants),len(times)))
        fitclock = np.zeros((len(ants),len(times)))
        fitoffset = np.zeros((len(ants),len(times)))
        fitrm = np.zeros((len(ants),len(times)))
        fitweights = np.ones((len(ants),len(times)))

        A = np.zeros((len(freqs), 2), dtype = float)
        A[:,0] = -8.44797245e9*2/freqs
        A[:,1] = freqs*2*np.pi*2
        #A[:,2] = np.ones((len(freqs),))
        sol = np.zeros((len(ants), 2))
        print ants
 
        fitguess = None
        fitrmguess = 0 # good guess
        for t, time in enumerate(times):

            # apply flags
            flags = (weights[0,:,:,t] == 0.) + (weights[1,:,:,t] == 0.)
            # phase.shape ~ (2, 36, 964, 4314) - pol, ant, freq, time
            phase       = phases[0,:,:,t] + phases[1,:,:,t] # not divide by 2 otherwise jump problem, then later fix this
            phase_diff  = phases[0,:,:,t] - phases[1,:,:,t] # not divide by 2 otherwise jump problem, then later fix this
            wav = c/freqs
            logging.debug('prep done')

            for refine in [1]:

                fitguess = initGuess(fitguess, ants, freqs, phase, refine)
                sol = fitting.fit(phase.T, A, fitguess, flags)

                logging.debug(str(sol.T))
                print A.shape, sol.shape
                logging.debug('done fit')

                # CHECK RESIDUALS
                residual = phase.T - np.dot(A, sol)
               # flags = 
                print "t:", t, "residual:", residual
                if residual < 0.5:
                    sol = figuess
                    break
               
                # high residual, do another cycle
#                logging.warning('Bad solution for ant: '+coord['ant']+' (time: '+str(t)+', resdiaul: '+str(residual)+', cycle: '+str(cycle)+')')
#                fitguess = None

                # Debug plot
                doplot = False
                if doplot:
                    if not 'matplotlib' in sys.modules:
                        import matplotlib as mpl
                        mpl.rc('font',size =8 )
                        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
                        mpl.use("Agg")
                    import matplotlib.pyplot as plt

                    fig = plt.figure()
                    fig.subplots_adjust(wspace=0)
                    ax = fig.add_subplot(110)

                    plotrm          = lambda RM, wav: np.mod( (2.*RM*wav*wav) + np.pi, 2.*np.pi)  -1.0*np.pi # notice the factor of 2
                    fitfuncfastplot = lambda p, freq: np.mod((4.*np.pi*p[1]*freq) - (2.*8.44797245e9*p[0]/freq) + (p[2]) + 1.0*np.pi, 2.*np.pi) - np.pi 

                    ax.plot(freq, np.mod(phase + np.pi, 2.*np.pi) - np.pi, 'or' )
                    ax.plot(freq, np.mod(phase_diff + np.pi, 2.*np.pi) - np.pi , '.', color='purple' )                           
    
                    TEC   = np.mod((-8.44797245e9*(2.*fitresult[0])/freq)+np.pi, 2*np.pi) - np.pi   # notice factor of 2 because rr+ll
                    Clock = np.mod((2.*np.pi*2.*fitresult[1]*freq )+np.pi, 2*np.pi) - np.pi   # notice factor of 2 because rr+ll
        
                    phase_total = (2.*np.pi*2.*fitresult[1]*freq)+(-8.44797245e9*(2.*fitresult[0])/freq)+fitresult[2]
                    residual    = np.mod(phase-phase_total+np.pi,2.*np.pi)-np.pi
                    ax.plot(freq, residual, '.', color='yellow')
        
                    idxl = int(min(freq)/1e4) 
                    idxh = int(max(freq)/1e4) 
                    bigfreqaxistmp = range(idxl, idxh)
                    bigfreqaxis    = np.array([float(i) for i in bigfreqaxistmp])
                    bigfreqaxis    = bigfreqaxis*1e4
            
                    ax.plot(bigfreqaxis, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")        
                    ax.plot(bigfreqaxis, plotrm(fitresultrm_wav, c/bigfreqaxis[:]), "-", color='purple')
            
                    ax.plot(freq, Clock, ',g') 
                    ax.plot(freq, TEC, ',b') 
                    ax.set_xlabel('freq')
                    ax.set_ylabel('phase')
                    ax.set_ylim(ymin=-np.pi, ymax=np.pi)
    
                    plt.savefig(str(t)+'_'+str(cycle)+'_'+coord['ant']+'.png', bbox_inches='tight')
          
            # rm is easy, go straight

            fittec[:,t] = sol[0]
            fitclock[:,t] = sol[1]
            fitoffset[:,t] = sol[2]
            fitrm[:,t] = solrm
            fitweights[:,t] = weight

       # return clock, tec, offset, rm
        return fittec, fitclock, fitoffset, fitrm, fitweights

def run( step, parset, H ):
    """
    Separate phase solutions into FR, Clock and TEC.

    The Clock and TEC values are stored in the specified output soltab with type 'clock', 'tec', 'FR'.
    """
    from losoto.h5parm import solFetcher, solWriter
    import numpy as np

    # get involved solsets using local step values or global values or all
    soltabs = getParSoltabs( step, parset, H )

    refAnt = parset.getString('.'.join(["LoSoTo.Steps", step, "RefAnt"]), '' )
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 1 )

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
            logging.error('Reference antenna '+refAnt+' not found.')
            return 1
        if refAnt == '': refAnt = ants[0]

        solsetname = soltabs[t].split('/')[0]
        st = H.makeSoltab(solsetname, 'tec',
                                 axesNames=['ant','time'], axesVals=[ants, times],
                                 vals=np.zeros((len(ants),len(times))),
                                 weights=np.ones((len(ants),len(times))))
        sw1 = solWriter(st)
        sw1.addHistory('Created by CLOCKTEC2 operation.')
        st = H.makeSoltab(solsetname, 'clock',
                                 axesNames=['ant','time'], axesVals=[ants, times],
                                 vals=np.zeros((len(ants),len(times))),
                                 weights=np.ones((len(ants),len(times))))
        sw2 = solWriter(st)
        sw2.addHistory('Created by CLOCKTEC2 operation.')
        st = H.makeSoltab(solsetname, 'rm',
                                 axesNames=['ant','time'], axesVals=[ants, times],
                                 vals=np.zeros((len(ants),len(times))),
                                 weights=np.ones((len(ants),len(times))))
        sw3 = solWriter(st)
        sw3.addHistory('Created by CLOCKTEC2 operation.')
            
        names = sf.getAxesNames()
        logging.debug('get vals start')
        vals, coord = sf.getValues(weight = False, retAxesVals = True)
        logging.debug('getweight start')
        weights = sf.getValues(weight = True, retAxesVals = False)

        if len(coord['freq']) < 10:
            logging.error('Clock/TEC separation needs at least 10 frequency channels, preferably distributed over a wide range.')
            return 1

        tec, clock, offset, rm, weight = fit_dTEC_dclock_dFR(np.squeeze(vals), np.squeeze(weights), coord)
    
        # remove pol from selection
        sw1.setSelection(ant=coord['ant'], time=coord['time'])
        sw2.setSelection(ant=coord['ant'], time=coord['time'])
        sw3.setSelection(ant=coord['ant'], time=coord['time'])
        sw1.setValues( tec )
        sw2.setValues( clock )
        sw3.setValues( rm )
        sw1.setValues( weight, weight=True )
        sw2.setValues( weight, weight=True )
        sw3.setValues( weight, weight=True )

    return 0
