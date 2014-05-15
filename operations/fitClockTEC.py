import numpy as np
import sys
import logging
from lofar.expion import baselinefitting  as fitting

def ClockTECfunc(xarray,par):
    delay=np.array([-1*par[1]*1e-9]).flatten() #in ns, array has dimension 1, even if scalar
    delayfact=2*np.pi*delay[:,np.newaxis]*xarray
    TEC=np.array([par[0]]).flatten();          # dTEC in TECU
    drefract=-8.4479745e9*TEC[:,np.newaxis]/xarray
    if len(par)>2:
        return drefract[:,np.newaxis,:]+delayfact+par[2]; #returns nTEC x nClock x nFreq
    
    return drefract[:,np.newaxis,:]+delayfact
    #return drefract+delayfact

def ClockTECfuncAllStations(xarray,par):
    delay=np.array([-1*par[1]*1e-9]).flatten() #in ns, array has dimension 1, even if scalar
    delayfact=2*np.pi*delay[:,np.newaxis]*xarray
    #print "delayfact",delayfact.shape
    TEC=np.array([par[0]]).flatten();          # dTEC in TECU
    drefract=-8.4479745e9*TEC[:,np.newaxis]/xarray
    #print "refract",drefract.shape
    if len(par)>2:
        return drefract+delayfact+par[2]; #returns nTEC x nClock x nFreq
    
    return drefract+delayfact
    #return drefract+delayfact

def getInitPar(data,dTECArray, dClockArray,freqs,ff=ClockTECfunc):
    '''initialize paramaters and unwraps data for fit'''
    if np.min(np.absolute(dTECArray))<1e-5:
        dTECArray+=0.0001 #just to prevent 0 to be there, since otherwise the fit might get stuck in 0 
    nT=dTECArray.shape[0]
    nD=dClockArray.shape[0]
    par=[dTECArray,dClockArray,0]
    # first check all unwrapping possibilities
    bigdata=ff(freqs,par) # returns array of shape nT,nD,nF
    wraps=np.around(np.divide(bigdata-data,2*np.pi));
    difference=  bigdata-data-wraps*2*np.pi
    offset=np.average(difference,axis=2);
    index=np.unravel_index(
        np.argmin(
            np.sum(np.absolute((difference.T-offset.T).T),axis=2))
        ,(nT,nD));
    OffsetIn=-1*offset[index];
    par=[dTECArray[index[0]],dClockArray[index[1]],OffsetIn];
    estimate=ff(freqs,par).flatten()
    wraps=np.around(np.divide(estimate-data,2*np.pi));
    data[:]=np.add(2*np.pi*wraps,data)
    return par

def getClockTECFit(ph,freq,stations,initSol=[],returnResiduals=True,chi2cut=1e8):
    nT=ph.shape[0]
    nF=freq.shape[0]
    nSt=ph.shape[2]
    data=ph
    tecarray=np.zeros((nT,nSt),dtype=np.float32)
    clockarray=np.zeros((nT,nSt),dtype=np.float32)
    if returnResiduals:
        residualarray=np.zeros((nT,nF,nSt),dtype=np.float32)
    nparms=2
    A=np.zeros((nF,nparms),dtype=np.float)
    A[:,1] = freq*2*np.pi*(-1e-9)
    A[:,0] = -8.44797245e9/freq
    base,steps=getPhaseWrapBase(freq)
    stepdTEC=np.abs(steps[0])*0.2
    stepDelay=np.abs(steps[1])*0.2
    succes=False
    initprevsol=False
    nrFail=0
    sol=np.zeros((nSt,nparms),dtype=np.float)
    for itm in range(nT):
        
        if itm%100==0 and itm>0:
            sys.stdout.write(str(itm)+'... '+str(sol[-1,0]-sol[0,0])+' '+str(sol[-1,1]-sol[0,1])+' '+str(sol[-1,-1]-sol[0,-1])+' ')
            sys.stdout.flush()

        if itm==0 or not succes:
         for ist in range(nSt):
             if itm==0 or not initprevsol:
                if hasattr(initSol,'__len__') and len(initSol)>ist:
                    iTEC1=initSol[ist,0]
                    iTEC2=initSol[ist,0]+stepdTEC
                    iD1=initSol[ist,1]
                    iD2=initSol[ist,1]+stepDelay
                else:
                 if 'CS' in stations[ist]:
                    iTEC1=-0.2
                    iTEC2=0.2
                    iD1=-4
                    iD2=4
                 else:
                    iTEC1=-1.5
                    iTEC2=1.5
                    iD1=-300
                    iD2=300
                logging.info("First %f %f %f %f %f %f "%(iTEC1,iTEC2,stepdTEC,iD1,iD2,stepDelay))

             else:
                
                iTEC1=prevsol[ist,0]-stepdTEC*nrFail
                iTEC2=prevsol[ist,0]+stepdTEC*(nrFail+1)
                if not fixedClockforCS or not 'CS' in stations[ist]: 
                    iD1=prevsol[ist,1]-stepDelay*nrFail
                    iD2=prevsol[ist,1]+stepDelay*(nrFail+1)
                else:
                    iD1=sol[ist,1]
                    iD2=sol[ist,1]+stepDelay
                    
                #print "Failure",iTEC1,iTEC2,iD1,iD2,nrFail

             dTECArray=np.arange(iTEC1,iTEC2,stepdTEC)
             dClockArray=np.arange(iD1,iD2,stepDelay)
             datatmp=ph[itm,:,ist]
             par = getInitPar(datatmp,dTECArray, dClockArray,freq,ClockTECfunc)
             sol[ist,:]=par[:nparms]
        for nr_iter in range(2):
            estimate=ClockTECfuncAllStations(freq,sol.T).reshape((nSt,nF)).T
            wraps=np.around(np.divide(estimate-data[itm],2*np.pi));
            data[itm,:]=np.add(2*np.pi*wraps,data[itm])
            sol=np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,data[itm])).T
        residual = data[itm] - np.dot(A, sol.T)
        residual = residual - residual[:, 0][:,np.newaxis]
        residual = np.remainder(residual+np.pi, 2*np.pi) - np.pi  
        if returnResiduals:
            residualarray[itm]=residual
        chi2=np.sum(np.square(np.degrees(residual)))/(nSt*nF)
        if chi2>chi2cut:
            logging.info("failure %f "%chi2 + str(sol))
            succes=False
            nrFail+=1
        else:
            prevsol=np.copy(sol)
#            print "succes",chi2,dTECArray.shape,dClockArray.shape
            succes=True
            initprevsol=True
            nrFail=0
        tecarray[itm]=sol[:,0]
        clockarray[itm]=sol[:,1]
    if returnResiduals:
        return tecarray,clockarray,residualarray
    return tecarray,clockarray

def getPhaseWrapBase(freqs):
    nF=freqs.shape[0]
    A=np.zeros((nF,2),dtype=np.float)
    A[:,1] = freqs*2*np.pi*(-1e-9)
    A[:,0] = -8.44797245e9/freqs
    steps=np.dot(np.dot(np.linalg.inv(np.dot(A.T,A)),A.T),2*np.pi*np.ones((nF,),dtype=np.float))
    basef=np.dot(A,steps)-2*np.pi
    return basef,steps

def getResidualPhaseWraps(avgResiduals,freqs):
    flags=avgResiduals[:,10]==0.
    nSt=avgResiduals.shape[1]
    nF=freqs.shape[0]
    wraps=np.zeros((nSt,),dtype=np.float)
    #tmpflags=np.sum(flags[:,np.sum(flags,axis=0)<(nF*0.5)],axis=1)
    tmpflags=flags
    tmpfreqs=freqs[np.logical_not(tmpflags)]
    tmpbasef,steps=getPhaseWrapBase(tmpfreqs)
    basef=np.zeros(freqs.shape)
    basef[np.logical_not(tmpflags)]=tmpbasef
    basef=basef.reshape((-1,1))
    
    data=avgResiduals[:,:]
        

    wraps=fitting.fit(data,basef,wraps,flags).flatten()
    return wraps,steps


def correctWraps(tecarray,residualarray,freq,pos):
    nT=tecarray.shape[0]
    nSt=tecarray.shape[1]
    avgResiduals=np.average(residualarray,axis=0)
    wraps,steps=getResidualPhaseWraps(avgResiduals,freq)
    lats=np.degrees(np.arctan2(pos[:,2],np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1])))
    lats-=lats[0]
    lons=np.degrees(np.arctan2(pos[:,1],pos[:,0]))
    lons-=lons[0]
    lonlat=np.concatenate((lons,lats)).reshape((2,)+lons.shape)
    for nr_iter in range(2):
        TEC=tecarray-tecarray[:,[0]]+steps[0]*(np.round(wraps)-np.round(wraps[0]))

        slope=np.dot(np.linalg.inv(np.dot(lonlat,lonlat.T)),np.dot(lonlat,TEC.T))
        chi2=np.sum(np.square(TEC-np.dot(lonlat.T,slope).T),axis=1)/nSt

        chi2select=chi2<np.average(chi2)
        chi2select=chi2<np.average(chi2[chi2select])
        offsets=-1*(np.average(TEC[chi2select]-np.dot(slope.T,lonlat)[chi2select],axis=0))*2.*np.pi/steps[0]
        remainingwraps=np.round(offsets/(2*np.pi))#-np.round(wraps[stationIndices])
        wraps+=remainingwraps

        if np.sum(np.absolute(remainingwraps))==0:
            continue

    return offsets,wraps,steps     


def doFit(phases,freqs,stations,station_positions,axes,refstIdx='superterp',stationSelect='BA',flagBadChannels=True,flagcut=1.5,chi2cut=30000.,removePhaseWraps=True,combine_pol=False,ignore_stations=["NOTHING_TO_IGNORE",]):
    #ake sure order of axes is as expected
    stidx=axes.index('ant')
    freqidx=axes.index('freq')
    timeidx=axes.index('time')
    polidx=axes.index('pol')
    data=phases.transpose((timeidx,freqidx,stidx,polidx))
    nT=data.shape[0]
    nF=data.shape[1]
    nSt=data.shape[2]
    npol=data.shape[3]
    if npol==4:
        data=data[:,:,:,(0,3)]
        npol=2

    if refstIdx=='superterp':
        superterpstations=[i for i in stations if i[:5] in ['CS002','CS003','CS004','CS005','CS006','CS007']]
        refstIdx=[i for (i,j) in enumerate(stations) if j in superterpstations ]
    if not hasattr(refstIdx,'__len__'):
        refstIdx=[refstIdx]

    #get phases from reference stations
    refdata=np.angle(np.average(np.cos(data[:,:,refstIdx,:])+1j*(np.sin(data[:,:,refstIdx,:])),axis=2))
    data-=refdata[:,:,np.newaxis]
    #flag bad channels
    if flagBadChannels:
        freqselect=np.ones((nF,),dtype=np.bool)
        for nr_iter in range(2):
            rms=np.std(np.std(np.std(data[:,:,:12],axis=0),axis=1),axis=1)
            freqselect=rms<flagcut*np.average(rms)
            logging.info("iter %d: flagging %d channels"%(nr_iter,np.sum(np.logical_not(freqselect))))
            freqs=freqs[freqselect]
            data=data[:,freqselect]
        nF=data.shape[1]
    #select stations
    if isinstance(stationSelect,str): 
        selectstations=[st for st in stations if stationSelect in st]   
    else:
        selectstations=list(stations[stationSelect])
    for ignore in ignore_stations:
        selectstations=[st for st in selectstations if not ignore in st]
    logging.info("%d selected stations: "%len(selectstations)+str(selectstations))
    stationIndices=np.array([idxst in selectstations for idxst in stations])
    data=data[:,:,stationIndices]
    stations=stations[stationIndices]
    nSt=data.shape[2]
    #cobine polarizationsif requested
    if combine_pol:
        if  npol==2:
            cdata=cos(data)+1j*sin(data)
            data=np.angle(np.average(cdata,axis=3)).reshape((nT,nF,nSt,1))
            npol=1
    #initialize arrays
    clock=np.zeros((nT,nSt,npol),dtype=np.float32)
    tec=np.zeros((nT,nSt,npol),dtype=np.float32)
    offset=np.zeros((nSt,npol),dtype=np.float32)
    for pol in range(npol):
        
        tecarray,clockarray,residualarray=getClockTECFit(data[:,:,:,pol],freqs,stations,initSol=[],returnResiduals=True,chi2cut=chi2cut)
        offset[:,pol],wraps,steps =correctWraps(tecarray,residualarray,freqs,station_positions)
        
        data[:,:,:,pol]+=offset[:,pol][np.newaxis,np.newaxis]
        initsol=np.zeros((nSt,2),dtype=np.float32)
        initsol[:,0]=tecarray[0,:]+wraps*steps[0]
        initsol[:,1]=clockarray[0,:]+wraps*steps[1]
        tecarray=0
        clockarray=0
        residualarray=0
        tec[:,:,pol],clock[:,:,pol]=getClockTECFit(data[:,:,:,pol],freqs,stations,initSol=initsol,returnResiduals=False)
        
    return clock,tec,offset,stations
