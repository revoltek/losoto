#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Clock/tec/FR separation module

import logging
from losoto.operations_lib import *
logging.debug('Loading CLOCKTEC2 module.')

c = 2.99792458e8

def fit_dFR(phases_rr, phases_ll, weights_rr, weights_ll, freq):


def fit_dTEC_dclock_dFR(phases_rr, phases_ll, weights_rr, weights_ll, freq, distance_station=1):
    freq_old = np.copy(freq)
    # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
    par3complex     = lambda p, freq, y: abs(np.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - np.cos(y)) + abs(np.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - np.sin(y))
    par2complex     = lambda p, freq, y: abs(np.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - np.cos(y)) + abs(np.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - np.sin(y))
    rmwavcomplex    = lambda RM, wav, y: abs(np.cos(2.*RM[0]*wav*wav)  - np.cos(y)) + abs(np.sin(2.*RM[0]*wav*wav)  - np.sin(y))

    # TODO: why?
    par3complex_w   = lambda p, freq, y: abs(np.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - np.cos(y))*(freq/1e5) + abs(np.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - np.sin(y))*(freq/1e5
    par2complex_w   = lambda p, freq, y: abs(np.cos((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - np.cos(y))*(freq/1e5) + abs(np.sin((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - np.sin(y))*(freq/1e5)
    rmwavcomplex_w  = lambda RM, wav, y: abs(np.cos(2.*RM[0]*wav*wav)  - np.cos(y))/wav + abs(np.sin(2.*RM[0]*wav*wav)  - np.sin(y))/wav

    # apply flags
    idx       = np.where(np.loagicanl_and(weights_rr != 0., weights_ll != 0.))
    freq      = freq[idx]
    phases_rr = phases_rr[idx]
    phases_ll = phases_ll[idx]
    logging.debug('Number of filtered out data points: '+str(len(idx[0]))

    if len(freq) >= 10: # prepare and make arrays if there is valid data 
        # TODO: why freq_old?
        freq      = (freq[0:len(freq_old)])
        phases_ll = (phases_ll[0: len(freq_old)])
        phases_rr = (phases_rr[0: len(freq_old)])
      
        phase       = (phases_rr + phases_ll)      # not divide by 2, then later fix this
        phase_diff  = (phases_rr - phases_ll)      # not divide by 2, then later fix this
      
        wav = c/freq
        pi = np.pi
        chi_old = 1.e9 
    else:
        logging.warning('No valid data found for clock/tec.')
        return (0.0, 0.0, 0.0, 0.0)
      
    # FIND INTIAL GUESS
    for dTEC in np.arange(-1.0,1.0, 0.01):
        for dclock in np.arange(-200e-9,200e-9,5e-9):
            phase_model = np.mod ( (4.*pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = np.mod (phase, 2*pi)
            angle       = pi - np.abs(np.abs(phase_model - phase_data) - pi)
            chi = np.sum(angle)

            if chi < chi_old:
                chi_old = chi
                fitguess     = [dclock,dTEC]
                #print 'Better fit', dclock, dTEC

    fitguess_1 = np.copy(fitguess)
    #print 'iter 1', fitguess

    for dTEC in np.arange(fitguess_1[1]-0.02,fitguess_1[1]+0.02, 0.002):
        for dclock in np.arange(fitguess_1[0]-8e-9,fitguess_1[0]+ 8e-9,1e-9):
            phase_model = np.mod ( (4.*pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
            phase_data  = np.mod (phase, 2*pi)
            angle       = pi - np.abs(np.abs(phase_model - phase_data) - pi)
            chi = np.sum(angle)

            if chi < chi_old:
                chi_old = chi
                fitguess = [dclock,dTEC]
                #print 'Better fit', dclock, dTEC
    #print 'iter 2', fitguess

    chi_old = 1e9
    for dFR in np.arange(-0.1,0.1,2e-4):
          phase_model = np.mod (2.*dFR*wav*wav, 2*pi) # notice the factor of 2
          phase_data  = np.mod (phase_diff, 2*pi)
          angle       = pi - np.abs(np.abs(phase_model - phase_data) - pi)
          chi         = np.sum(angle)

          if chi < chi_old:
              chi_old = chi
              fitrmguess = dFR
              #print 'Better fit', fitrmguess

    fitrmguess_1 = np.copy(fitrmguess)
    for dFR in np.arange(fitrmguess_1-5e-4,fitrmguess_1+5e-4,0.5e-5):
        phase_model = np.mod (2.*dFR*wav*wav, 2*pi) # notice the factor of 2
        phase_data  = np.mod (phase_diff, 2*pi)
        angle       = pi - np.abs(np.abs(phase_model - phase_data) - pi)
        chi         = np.sum(angle)

        if chi < chi_old:
            chi_old = chi
            fitrmguess = dFR
            #print 'Better fit', fitrmguess

    # DO THE FITTING 
    # SOLVE Clock-TEC anticorrelation problem on short baselines             
    freq = freq.astype(np.float64)
    phase = phase.astype(np.float64) #epsfcn=1e-7
      
    # TODO: what is distant station?
    if distance_station < 0. :   #15.0*1e3: DOES NOT WORK, NEED 3 par FIT
        fitresult, success  = scipy.optimize.leastsq(par2complex, fitguess, args=(freq, phase))
        #fitresult = fitguess
    else:
        fitresult, success   = scipy.optimize.leastsq(par3complex, [fitguess[0], fitguess[1], 0.0], args=(freq, phase), maxfev=10000)
        #fitresult = [fitguess[0], fitguess[1], 0.0]
        #print fitresult, success
        fitresultrm_wav, success = scipy.optimize.leastsq(rmwavcomplex, [fitrmguess], args=(wav, phase_diff))
              
    show_plot = False
    if show_plot:  

        plotrm          = lambda RM, wav: np.mod( (2.*RM*wav*wav) +1.0*pi, 2.*pi)  -1.0*pi # notice the factor of 2
        fitfuncfastplot = lambda p, freq: np.mod((4.*pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + (p[2])+ 1.0*pi, 2.*pi) -1.0*pi 
        
        #if len(fitresult ==2): 
        #  fitresult =[fitresult[0], fitresult[1],0]
        #  print 'Here'
        #  #fitresult = fitresult #   [fitguess[0],  fitguess[1], 0]
        
        #fitresult = [[fitresultguess[0],fitresultguess[0],0.0]]
        
        matplotlib.pyplot.plot(freq, np.mod(phase + 1.0*pi, 2.*pi) -1.0*pi, 'or' )
        matplotlib.pyplot.plot(freq, np.mod(phase_diff + 1.0*pi, 2.*pi) -1.0*pi , '.', color='purple' )                           

        TEC   = np.mod((-8.44797245e9*(2.*fitresult[1])/freq)+np.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
        Clock = np.mod((2.*np.pi*   2.*fitresult[0]*freq )+np.pi, 2*pi) - pi   # notice factor of 2 because rr+ll
    
        phase_total = (2.*np.pi*2.*fitresult[0]*freq)+(-8.44797245e9*(2.*fitresult[1])/freq)+fitresult[2]
        residual    = np.mod(phase-phase_total+pi,2.*pi)-pi
        matplotlib.pyplot.plot(freq, residual, '.', color='yellow')
    
        idxl = int(min(freq_old)/1e4) 
        idxh = int(max(freq_old)/1e4) 
        bigfreqaxistmp = range(idxl, idxh)
        bigfreqaxis    =  np.array([float(i) for i in bigfreqaxistmp])
        bigfreqaxis    = bigfreqaxis*1e4
        
        matplotlib.pyplot.plot (bigfreqaxis, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")        
        matplotlib.pyplot.plot (bigfreqaxis, plotrm(fitresultrm_wav, c/bigfreqaxis[:]), "-", color='purple')
        #matplotlib.pyplot.plot (freq, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")        
        
        matplotlib.pyplot.plot(freq, Clock, ',g') 
        matplotlib.pyplot.plot(freq, TEC, ',b') 
        matplotlib.pyplot.xlabel('freq')
        matplotlib.pyplot.ylabel('phase')

        matplotlib.pyplot.show()
      
    # return clock, tec, offset, rm
    return (fitresult[0], fitresult[1], fitresult[2], fitresultrm_wav)

def run( step, parset, H ):
    """
    Separate phase solutions into FR, Clock and TEC.

    The Clock and TEC values are stored in the specified output soltab with type 'clock', 'tec', 'FR'.
    """
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter

    # get involved solsets using local step values or global values or all
    soltabs = getParSoltabs( step, parset, H )

    flagBadChannels = parset.getBool('.'.join(["LoSoTo.Steps", step, "FlagBadChannels"]), True )
    flagCut = parset.getFloat('.'.join(["LoSoTo.Steps", step, "FlagCut"]), 1.5 )
    chi2cut = parset.getFloat('.'.join(["LoSoTo.Steps", step, "Chi2cut"]), 30000. )
    refAnt = parset.getFloat('.'.join(["LoSoTo.Steps", step, "RefAnt"]), '' )

    for soltab in openSoltabs( H, soltabs ):
        logging.info("--> Working on soltab: "+soltab._v_name)
        sf = solFetcher(soltab)

        # this will make a selection for the getValues() and getValuesIter()
        userSel = {}
        for axis in t.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        times = sf.getAxisValues('time')
        ants = sf.getAxisValues('ant')

        # some checks
        solType = sf.getType()
        if solType != 'phase':
           logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase. Ignoring.")
           continue

        if refAnt != '' and not refAnt in ants:
            logging.error('Reference antenna '+refAnt+' not found.')
            return 1
        if refAnt == '': refAnt = ants[0]

        st = H.makeSoltab(solsetname, 'clock',
                                 axesNames=['time', 'ant'], axesVals=[times, ants],
                                 vals=np.zeros((len(times),len(ants))),
                                 weights=np.ones((len(times),len(ants))))
        sw1 = solWriter(st)
        sw1.addHistory('Created by CLOCKTEC2 operation.')
        st = H.makeSoltab(solsetname, 'tec',
                                 axesNames=['time', 'ant'], axesVals=[times, ants],
                                 vals=np.zeros((len(times),len(ants))),
                                 weights=np.ones((len(times),len(ants))))
        sw2 = solWriter(st)
        sw2.addHistory('Created by CLOCKTEC2 operation.')
        st = H.makeSoltab(solsetname, 'rm',
                                 axesNames=['time','ant'], axesVals=[times, ants],
                                 vals=np.zeros((len(times),len(ants))),
                                 weights=np.ones((len(times),len(ants))))
        sw3 = solWriter(st)
        sw3.addHistory('Created by CLOCKTEC2 operation.')
            
        import time
        start = time.time()
        for vals, weights, coord, selection in t.getValuesIter(returnAxes=['ant','freq','pol'], weight=True, reference = refAnt):
            print 'read', start - time.time()
            start = time.time()

            if len(coord['ant']) < 2:
                logging.error('Clock/TEC separation needs at least 2 antennas selected.')
                return 1
            if len(coord['freq']) < 10:
                logging.error('Clock/TEC separation needs at least 10 frequency channels, preferably distributed over a wide range.')
                return 1
            
            print vals.shape
            sys.exit(1)

            clock = np.zeros_like(vals)
            tec = np.zeros_like(vals)
            offset = np.zeros_like(vals)
            rm = np.zeros_like(vals)
            print 'allocate', start - time.time()
            start = time.time()

            for i, ant in enumerate(coord['ant']):
                if ant == refAnt: continue
                clock[i:], tec[i:], offset[i:], rm[i:] = fit_dTEC_dclock_dFR(vals[i:], vals[i:], weights[i:], weights[i:], coord['freq'])

            print 'ct sep', start - time.time()
            start = time.time()
            
            sw1.selection = selection
            sw2.selection = selection
            sw3.selection = selection
            sw1.setValues(clock)
            sw2.setValues(tec)
            sw3.setValues(rm)
            
            print 'write', start - time.time()
            start = time.time()

    return 0
