#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Clock/tec/FR separation module

import logging
from losoto.operations_lib import *
logging.debug('Loading CLOCKTEC2 module.')

import multiprocessing
inQueue = multiprocessing.JoinableQueue()
outQueue = multiprocessing.Queue()

c = 2.99792458e8

#def fit_dFR(phases_rr, phases_ll, weights_rr, weights_ll, freq):

class multiThread(multiprocessing.Process):
    """
    This class is a working thread which load parameters from a queue and
    run the flagging on a chunk of data
    """

    def __init__(self, inQueue, outQueue):
        multiprocessing.Process.__init__(self)
        self.inQueue = inQueue
        self.outQueue = outQueue

    def run(self):

        while True:
            parms = self.inQueue.get()

            # poison pill
            if parms is None:
                self.inQueue.task_done()
                break

            self.fit_dTEC_dclock_dFR(*parms)
            self.inQueue.task_done()

    def fit_dTEC_dclock_dFR(self, phases_rr, phases_ll, weights_rr, weights_ll, freq, ant):
        import numpy as np
        import scipy.optimize
    
        # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
        par3complex = lambda p, freq, y: abs(np.cos((4.*np.pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - np.cos(y)) + abs(np.sin((4.*np.pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + p[2]) - np.sin(y))
    #    par2complex = lambda p, freq, y: abs(np.cos((4.*np.pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - np.cos(y)) + abs(np.sin((4.*np.pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq)) - np.sin(y))
        rmwavcomplex = lambda RM, wav, y: abs(np.cos(2.*RM[0]*wav*wav)  - np.cos(y)) + abs(np.sin(2.*RM[0]*wav*wav)  - np.sin(y))
    
        # apply flags
        idx       = np.where(np.logical_and(weights_rr != 0., weights_ll != 0.))
        freq      = freq[idx]
        phases_rr = phases_rr[idx]
        phases_ll = phases_ll[idx]
        #logging.debug('Number of filtered out data points: '+str(len(weights_rr) - len(idx[0])))
    
        if len(freq) < 10:
            logging.warning('No valid data found for clock/tec.')
            return (0.0, 0.0, 0.0, 0.0)
    
        phase       = (phases_rr + phases_ll)      # not divide by 2, then later fix this
        phase_diff  = (phases_rr - phases_ll)      # not divide by 2, then later fix this
        wav = c/freq

        # FIND INTIAL GUESS
        chis = []; dTECs = []; dclocks = []
        for dTEC in np.arange(-0.8,0.8, 0.01):
            for dclock in np.arange(-100e-9,100e-9,5e-9):
                phase_model = np.mod( (4.*np.pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*np.pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
                phase_data  = np.mod(phase, 2*np.pi)
                angle       = np.pi - np.abs(np.abs(phase_model - phase_data) - np.pi)
                chis.append(np.sum(angle))
                dTECs.append(dTEC)
                dclocks.append(dclock)
    
        idx = chis.index(min(chis))
        fitguess = [dTECs[idx], dclocks[idx]]
    
        chis = []; dTECs = []; dclocks = []
        for dTEC in np.arange(fitguess[0]-0.02, fitguess[0]+0.02, 0.002):
            for dclock in np.arange(fitguess[1]-8e-9, fitguess[1]+ 8e-9, 1e-9):
                phase_model = np.mod ( (4.*np.pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*np.pi)  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
                phase_data  = np.mod (phase, 2*np.pi)
                angle       = np.pi - np.abs(np.abs(phase_model - phase_data) - np.pi)
                chis.append(np.sum(angle))
                dTECs.append(dTEC)
                dclocks.append(dclock)
    
        idx = chis.index(min(chis))
        fitguess = [dTECs[idx], dclocks[idx]]
    
        chis = []
        for dFR in np.arange(-0.1, 0.1, 2e-4):
              phase_model = np.mod (2.*dFR*wav*wav, 2*np.pi) # notice the factor of 2
              phase_data  = np.mod (phase_diff, 2*np.pi)
              angle       = np.pi - np.abs(np.abs(phase_model - phase_data) - np.pi)
              chis.append(np.sum(angle))

        idx = chis.index(min(chis))
        fitrmguess = np.arange(-0.1, 0.1, 2e-4)[idx]
    
        chis = []
        for dFR in np.arange(fitrmguess-5e-4,fitrmguess+5e-4,0.5e-5):
            phase_model = np.mod (2.*dFR*wav*wav, 2*np.pi) # notice the factor of 2
            phase_data  = np.mod (phase_diff, 2*np.pi)
            angle       = np.pi - np.abs(np.abs(phase_model - phase_data) - np.pi)
            chis.append(np.sum(angle))

        idx = chis.index(min(chis))
        fitrmguess = np.arange(fitrmguess-5e-4,fitrmguess+5e-4,0.5e-5)[idx]
    
        # DO THE FITTING 
        # SOLVE Clock-TEC anticorrelation problem on short baselines             
        freq = freq.astype(np.float64)
        phase = phase.astype(np.float64) #epsfcn=1e-7
          
        fitresult, success = scipy.optimize.leastsq(par3complex, [fitguess[0], fitguess[1], 0.0], args=(freq, phase), maxfev=10000)
        fitresultrm_wav, success = scipy.optimize.leastsq(rmwavcomplex, [fitrmguess], args=(wav, phase_diff))
                  
        show_plot = False
        if show_plot:  
    
            plotrm          = lambda RM, wav: np.mod( (2.*RM*wav*wav) +1.0*np.pi, 2.*np.pi)  -1.0*np.pi # notice the factor of 2
            fitfuncfastplot = lambda p, freq: np.mod((4.*np.pi*p[0]*freq) - (2.*8.44797245e9*p[1]/freq) + (p[2])+ 1.0*np.pi, 2.*np.pi) -1.0*np.pi 
            
            matplotlib.pyplot.plot(freq, np.mod(phase + 1.0*np.pi, 2.*np.pi) -1.0*np.pi, 'or' )
            matplotlib.pyplot.plot(freq, np.mod(phase_diff + 1.0*np.pi, 2.*np.pi) -1.0*np.pi , '.', color='purple' )                           
    
            TEC   = np.mod((-8.44797245e9*(2.*fitresult[1])/freq)+np.pi, 2*np.pi) - np.pi   # notice factor of 2 because rr+ll
            Clock = np.mod((2.*np.pi*   2.*fitresult[0]*freq )+np.pi, 2*np.pi) - np.pi   # notice factor of 2 because rr+ll
        
            phase_total = (2.*np.pi*2.*fitresult[0]*freq)+(-8.44797245e9*(2.*fitresult[1])/freq)+fitresult[2]
            residual    = np.mod(phase-phase_total+np.pi,2.*pnp.i)-np.pi
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
        #return (fitresult[0], fitresult[1], fitresult[2], fitresultrm_wav)
        self.outQueue.put([fitresult[0], fitresult[1], fitresult[2], fitresultrm_wav[0], ant])

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

    # start processes for multi-thread
    logging.debug('Spowning %i threads...' % ncpu)
#    for i in xrange(ncpu):
#        t = multiThread(inQueue, outQueue)
#        t.start()

    for t, soltab in enumerate(openSoltabs( H, soltabs )):
        logging.info("--> Working on soltab: "+soltab._v_name)
        sf = solFetcher(soltab)

        # this will make a selection for the getValues() and getValuesIter()
        userSel = {}
        for axis in sf.getAxesNames():
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

        solsetname = soltabs[t].split('/')[0]
        st = H.makeSoltab(solsetname, 'clock',
                                 axesNames=['ant','time'], axesVals=[ants, times],
                                 vals=np.zeros((len(ants),len(times))),
                                 weights=np.ones((len(ants),len(times))))
        sw1 = solWriter(st)
        sw1.addHistory('Created by CLOCKTEC2 operation.')
        st = H.makeSoltab(solsetname, 'tec',
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
            
        import time
        start = time.time()
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=['ant','freq','pol'], weight=True, reference = refAnt):
            print 'read', time.time() - start
            start = time.time()

            if len(coord['freq']) < 10:
                logging.error('Clock/TEC separation needs at least 10 frequency channels, preferably distributed over a wide range.')
                return 1
            
            # assume standard ['pol', 'dir', 'ant', 'freq', 'time'] order
            names = sf.getAxesNames()
            assert names == ['pol', 'dir', 'ant', 'freq', 'time']
        
            clock  = np.zeros((len(coord['ant'])))
            tec    = np.zeros((len(coord['ant'])))
            offset = np.zeros((len(coord['ant'])))
            rm     = np.zeros((len(coord['ant'])))

            for i in xrange(ncpu):
                t = multiThread(inQueue, outQueue)
                t.start()

            for i, ant in enumerate(coord['ant']):
                if ant == refAnt: continue
                inQueue.put([vals[0,i,:], vals[1,i,:], weights[0,i,:], weights[1,i,:], coord['freq'], i])
                #clock[i], tec[i], offset[i], rm[i] = fit_dTEC_dclock_dFR(vals[0,i,:], vals[1,i,:], weights[0,i,:], weights[1,i,:], coord['freq'])

            # add poison pills to kill processes
            for i in xrange(ncpu):
                inQueue.put(None)

            print 'prep', time.time() - start
            start = time.time()
            # wait for all jobs to finish
            inQueue.join()

            print 'ct sep', time.time() - start
            start = time.time()

            for i in xrange(len(coord['ant'])-1):
                q = outQueue.get()
                ant = q[-1]
                clock[ant], tec[ant], offset[ant], rm[ant], ant = q
    
            # remove pol from selection
            sw1.selection = [ selection[names.index('ant')], selection[names.index('time')] ]
            sw2.selection = [ selection[names.index('ant')], selection[names.index('time')] ]
            sw3.selection = [ selection[names.index('ant')], selection[names.index('time')] ]
            sw1.setValues( np.expand_dims(clock, axis=1) )
            sw2.setValues( np.expand_dims(tec, axis=1) )
            sw3.setValues( np.expand_dims(rm, axis=1) )
            
            print 'write', time.time() - start
            start = time.time()

    return 0
