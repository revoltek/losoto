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


    def initGuessMin(self, dTECrange, dclockrange, phase):
        pass
        

    def initGuess(self, fitguess, ant, freq, phase, refine=1):
        import numpy as np

        # FIND INTIAL GUESS
        if fitguess == None:
            chis = []; dTECs = []; dclocks = []
            if 'CS' in ant:
                dclockrange = np.arange(-10e-9, 10e-9, 1.e-9/refine)
                dTECrange = np.arange(-0.01, 0.01, 0.001/refine)
            else: 
                dclockrange = np.arange(-200e-9, 200e-9, 5.e-9/refine)
                dTECrange = np.arange(-0.5, 0.5, 0.01/refine)

            for dTEC in dTECrange:
                for dclock in dclockrange:
                    phase_model = np.mod( (4.*np.pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*np.pi )  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
                    phase_data  = np.mod( phase, 2*np.pi )
                    angle       = np.pi - np.abs(np.abs(phase_model - phase_data) - np.pi)
                    chis.append(np.sum(angle))
                    dTECs.append(dTEC)
                    dclocks.append(dclock)
    
            idx = chis.index(min(chis))
            fitguess = [dTECs[idx], dclocks[idx]]
    
        chis = []; dTECs = []; dclocks = []
        if 'CS' in ant:
            dclockrange = np.arange(fitguess[1]-1.5e-9, fitguess[1]+1.5e-9, 0.1e-9/refine)
            dTECrange = np.arange(fitguess[0]-0.0015, fitguess[0]+0.0015, 0.0001/refine)
        else: 
            dclockrange = np.arange(fitguess[1]-7.5e-9, fitguess[1]+7.5e-9, 1.e-9/refine)
            dTECrange = np.arange(fitguess[0]-0.015, fitguess[0]+0.015, 0.001/refine)

        for dTEC in dTECrange:
            for dclock in dclockrange:
                phase_model = np.mod( (4.*np.pi*dclock*freq) - (2.*8.44797245e9*dTEC/freq), 2*np.pi )  # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
                phase_data  = np.mod( phase, 2*np.pi )
                angle       = np.pi - np.abs(np.abs(phase_model - phase_data) - np.pi)
                chis.append(np.sum(angle))
                dTECs.append(dTEC)
                dclocks.append(dclock)
    
        idx = chis.index(min(chis))
        return [dTECs[idx], dclocks[idx]]
 

    def fit_dTEC_dclock(self, phases, weights, coord):
        import numpy as np
        import scipy.optimize
    
        # NOTE THE *2 to use rr+ll instead of 0.5*(rr+ll)
        par3complex = lambda p, freq, y: abs(np.cos((4.*np.pi*p[1]*freq) - (2.*8.44797245e9*p[0]/freq)) - np.cos(y)) + abs(np.sin((4.*np.pi*p[1]*freq) - (2.*8.44797245e9*p[0]/freq)) - np.sin(y))
#        par3complex = lambda p, freq, y: abs(np.cos((4.*np.pi*p[1]*freq) - (2.*8.44797245e9*p[0]/freq) + p[2]) - np.cos(y)) + abs(np.sin((4.*np.pi*p[1]*freq) - (2.*8.44797245e9*p[0]/freq) + p[2]) - np.sin(y))

        times = np.copy(coord['time'])
        fittec = np.zeros(len(times))
        fitclock = np.zeros(len(times))
        fitoffset = np.zeros(len(times))
        fitrm = np.zeros(len(times))
        fitweights = np.ones(len(times))
    
        fitguess = None
        for t, time in enumerate(times):

            # apply flags
            idx       = np.where(np.logical_and(weights[0,:,t] != 0., weights[1,:,t] != 0.))
            freq      = np.copy(coord['freq'])[idx]
            phase_rr = phases[0,:,t][idx]
            phase_ll = phases[1,:,t][idx]
            #logging.debug('Number of filtered out data points: '+str(len(weights_rr) - len(idx[0])))
        
            if len(freq) < 10:
                logging.warning('No valid data found for clock/tec.')
                continue
        
            phase       = (phase_rr + phase_ll)      # not divide by 2 otherwise jump problem, then later fix this
    
            weight = 0
            cycle = 0
            for refine in [1]:#[1,1,2,3,4]:

                fitguess = self.initGuess(fitguess, coord['ant'], freq, phase, refine)
        
                # DO THE FITTING 
                freq = freq.astype(np.float64)
                phase = phase.astype(np.float64)
              
                #fitresult, success = scipy.optimize.leastsq(par3complex, [fitguess[0], fitguess[1], 0.0], args=(freq, phase), maxfev=10000)
                fitresult, success = scipy.optimize.leastsq(par3complex, [fitguess[0], fitguess[1]], args=(freq, phase), maxfev=10000)
    
                # CHECK RESIDUALS
                phase_total = (2.*np.pi*2.*fitresult[1]*freq) + (-8.44797245e9*(2.*fitresult[0])/freq)
                #phase_total = (2.*np.pi*2.*fitresult[1]*freq) + (-8.44797245e9*(2.*fitresult[0])/freq) + fitresult[2]
                residual = np.mean(abs(np.mod(phase-phase_total+np.pi,2.*np.pi)-np.pi))

                print coord['ant'], fitresult
                #print 'guess', fitguess, fitrmguess
                #print "t:", t, "cycle:", cycle, "residual:", residual

                if residual < 0.5:
                    fitguess = fitresult
                    weight = 1
                    break
               
                # high residual, do another cycle
                logging.warning('Bad solution for ant: '+coord['ant']+' (time: '+str(t)+', resdiaul: '+str(residual)+', cycle: '+str(cycle)+')')
                fitguess = None

                # Debug plot
                doplot = False
                if doplot:
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

                    fitfuncfastplot = lambda p, freq: np.mod((4.*np.pi*p[1]*freq) - (2.*8.44797245e9*p[0]/freq) + 1.0*np.pi, 2.*np.pi) - np.pi 

                    ax.plot(freq, np.mod(phase + np.pi, 2.*np.pi) - np.pi, 'or' )
                    ax.plot(freq, np.mod(phase_rr + np.pi, 2.*np.pi) - np.pi, '.b' )
                    ax.plot(freq, np.mod(phase_ll + np.pi, 2.*np.pi) - np.pi, '.g' )
    
                    TEC   = np.mod((-8.44797245e9*(2.*fitresult[0])/freq)+np.pi, 2*np.pi) - np.pi   # notice factor of 2 because rr+ll
                    Clock = np.mod((2.*np.pi*2.*fitresult[1]*freq )+np.pi, 2*np.pi) - np.pi   # notice factor of 2 because rr+ll
        
                    phase_total = (2.*np.pi*2.*fitresult[1]*freq)+(-8.44797245e9*(2.*fitresult[0])/freq)#+fitresult[2]
                    residual    = np.mod(phase-phase_total+np.pi,2.*np.pi)-np.pi
                    ax.plot(freq, residual, '.', color='yellow')
        
                    idxl = int(min(freq)/1e4) 
                    idxh = int(max(freq)/1e4) 
                    bigfreqaxistmp = range(idxl, idxh)
                    bigfreqaxis    = np.array([float(i) for i in bigfreqaxistmp])
                    bigfreqaxis    = bigfreqaxis*1e4
            
                    ax.plot(bigfreqaxis, fitfuncfastplot(fitresult, bigfreqaxis[:]), "r-")        
            
                    ax.plot(freq, Clock, ',g') 
                    ax.plot(freq, TEC, ',b') 
                    ax.set_xlabel('freq')
                    ax.set_ylabel('phase')
                    ax.set_ylim(ymin=-np.pi, ymax=np.pi)
    
                    logging.warning('Save pic: '+str(t)+'_'+str(cycle)+'_'+coord['ant']+'.png')
                    plt.savefig(str(t)+'_'+str(cycle)+'_'+coord['ant']+'.png', bbox_inches='tight')
                    del fig

                cycle += 1
          
            fittec[t] = fitresult[0]
            fitclock[t] = fitresult[1]
            fitoffset[t] = 0 #fitresult[2]
            fitweights[t] = weight

       # return clock, tec, offset, rm
        self.outQueue.put([fittec, fitclock, fitoffset, fitweights, coord])

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
    for i in xrange(ncpu):
        t = multiThread(inQueue, outQueue)
        t.start()

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
            
        import time
        proc=0
        for vals, weights, coord, selection in sf.getValuesIter(returnAxes=['freq','pol','time'], weight=True, reference = refAnt):

            if len(coord['freq']) < 10:
                logging.error('Clock/TEC separation needs at least 10 frequency channels, preferably distributed over a wide range.')
                return 1

            if not coord['ant'] == refAnt:
                inQueue.put([vals, weights, coord])
                proc += 1

        # TODO: add a check on axis order

        # add poison pills to kill processes
        for i in xrange(ncpu):
            inQueue.put(None)

        start = time.time()
        # wait for all jobs to finish
        inQueue.join()
        print 'ct sep', time.time() - start
        start = time.time()

        names = sf.getAxesNames()
        for i in xrange(proc):
            q = outQueue.get()
            tec, clock, offset, weight, coord = q
    
            # remove pol from selection
            sw1.setSelection(ant=coord['ant'], time=coord['time'])
            sw2.setSelection(ant=coord['ant'], time=coord['time'])
            sw1.setValues( np.expand_dims(tec, axis=1) )
            sw2.setValues( np.expand_dims(clock, axis=1) )
            sw1.setValues( np.expand_dims(weight, axis=1), weight=True )
            sw2.setValues( np.expand_dims(weight, axis=1), weight=True )

    return 0
