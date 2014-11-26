import numpy as np
import numpy.ma as ma
import sys
import logging
from lofar.expion import baselinefitting  as fitting

def ClockTECfunc(xarray,par):
    delay=np.array([-1*par[1]*1e-9]).flatten() #in ns, array has dimension 1, even if scalar
    delayfact=2*np.pi*delay[:,np.newaxis]*xarray
    TEC=np.array([par[0]]).flatten();          # dTEC in TECU
    drefract=-8.4479745e9*TEC[:,np.newaxis]/xarray
    if len(par)>2:
        return drefract[:,np.newaxis,:]+delayfact[np.newaxis]+par[2] #returns nTEC x nClock x nFreq
    
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
        return drefract+delayfact+par[2][:,np.newaxis] #returns nTEC x nClock x nFreq
    
    return drefract+delayfact
    #return drefract+delayfact



def getInitClock(data,freq):
    nF=freq.shape[0]
    avgdata=np.ma.sum(np.cos(data)+1.j*np.sin(data),axis=0).swapaxes(0,-2)
    avgdata=np.arctan2(np.imag(avgdata),np.real(avgdata))
    nSt=avgdata.shape[0]
    npol=avgdata.shape[2]
    for ist in range(nSt):
        for pol in range(npol):
            mymask=avgdata[ist,:,pol].mask
            if not hasattr(mymask,'__len__'):
                mymask=np.ones(avgdata[ist,:,pol].shape,dtype=bool)*mymask
            if avgdata[ist,:,pol].count()<1:
                avgdata[ist,:,pol].mask[0]=False
            #logging.info("mask station %d pol %d "%(ist,pol) +str(mymask)) 
            #logging.info("average data station %d pol %d "%(ist,pol) +str(avgdata[ist,:,pol])) 
            avgdata[ist,:,pol][~mymask]=np.float32(np.unwrap(avgdata[ist,:,pol][~mymask]))
            #logging.info("average unwrapped data station %d pol %d "%(ist,pol) +str(avgdata[ist,:,pol])) 
            #logging.info("remainder " +str(np.remainder(avgdata[ist,:,pol]+np.pi,2*np.pi)-np.pi)) 
    #print data.shape,avgdata.shape #stationsx freq x pol
    A=np.ones((nF,2),dtype=np.float)
    A[:,1] = freq*2*np.pi*(-1e-9)
    return np.ma.dot(np.linalg.inv(np.dot(A.T,A)),np.ma.dot(A.T,avgdata).swapaxes(0,-2))

def getInitPar(data,dTECArray, dClockArray,freqs,ff=ClockTECfunc):
    '''initialize paramaters and unwraps data for fit'''
    if np.min(np.absolute(dTECArray))<1e-5:
        dTECArray+=0.0001 #just to prevent 0 to be there, since otherwise the fit might get stuck in 0 
    nT=dTECArray.shape[0]
    nD=dClockArray.shape[0]
    par=[dTECArray,dClockArray,0]
    # first check all unwrapping possibilities
    bigdata=ff(freqs,par) # returns array of shape nT,nD,nF
    wraps=np.ma.around(np.divide(bigdata-data[np.newaxis,np.newaxis],2*np.pi));
    difference=  bigdata-data[np.newaxis,np.newaxis]-wraps*2*np.pi
    index=np.unravel_index(
        np.ma.argmin(
            np.ma.var(difference,axis=2))
        ,(nT,nD))
    OffsetIn=-1*np.ma.mean(difference[index]);
    par=[dTECArray[index[0]],dClockArray[index[1]],OffsetIn];
    estimate=ff(freqs,par).flatten()
    #logging.info("estimate "+str(estimate))
    wraps=np.ma.around(np.divide(estimate-data,2*np.pi));
    data[:]=np.add(2*np.pi*wraps,data)
    return par

def getClockTECFit(ph,freq,stations,initSol=[],returnResiduals=True,chi2cut=1e8 ,nparms=2):
    stepFraction=0.2
    nT=ph.shape[0]
    nF=freq.shape[0]
    nSt=ph.shape[2]
    data=ph
    logging.info("fitting masked data "+str(ph.count(axis=0)))
    tecarray=np.zeros((nT,nSt),dtype=np.float32)
    clockarray=np.zeros((nT,nSt),dtype=np.float32)
    if nparms>2:
        nparms=3
        fitoffsetarray=np.zeros((nT,nSt),dtype=np.float32)
    if returnResiduals:
        residualarray=np.zeros((nT,nF,nSt),dtype=np.float32)
    A=np.ones((nF,nparms),dtype=np.float)
    A[:,1] = freq*2*np.pi*(-1e-9)
    A[:,0] = -8.44797245e9/freq
    base,steps=getPhaseWrapBase(freq)
    stepdTEC=np.abs(steps[0])*stepFraction
    stepDelay=np.abs(steps[1])*stepFraction
    succes=False
    initprevsol=np.zeros(nSt,dtype=bool)
    nrFail=np.zeros(nSt,dtype=int)
    sol=np.zeros((nSt,nparms),dtype=np.float)
    prevsol=np.copy(sol)
    for itm in range(nT):
        
        if itm%100==0 and itm>0:
            sys.stdout.write(str(itm)+'... '+str(sol[-1,0]-sol[0,0])+' '+str(sol[-1,1]-sol[0,1])+' '+str(sol[-1,-1]-sol[0,-1])+' ')
            sys.stdout.flush()

        if itm==0 or not succes:
         for ist in range(nSt):
             if itm==0 or not initprevsol[ist]:
                if hasattr(initSol,'__len__') and len(initSol)>ist:
                    iTEC1=initSol[ist,0]
                    iTEC2=initSol[ist,0]+stepdTEC
                    iD1=initSol[ist,1]
                    iD2=initSol[ist,1]+stepDelay
                else:
                 if 'CS' in stations[ist]:
                    iTEC1=-0.2
                    iTEC2=0.2
                    iD1=-20
                    iD2=20
                 else:
                     if 'RS' in stations[ist]:
                         iD1=-150
                         iD2=150
                         iTEC1=-1.5
                         iTEC2=1.5
                     else: #large TEC variation for EU stations
                         iD1=-250
                         iD2=250
                         iTEC1=-5
                         iTEC2=5
                        
                logging.info("First %f %f %f %f %f %f "%(iTEC1,iTEC2,stepdTEC,iD1,iD2,stepDelay))

             else:
                sol[ist,:]=prevsol[ist,:]
                iTEC1=prevsol[ist,0]-min(1.5,stepdTEC*int(nrFail[ist]/1))#/stepFraction
                iTEC2=prevsol[ist,0]+min(1.5,stepdTEC*(int(nrFail[ist]/1)+1))#/stepFraction
                if not 'CS' in stations[ist]: 
                    iD1=prevsol[ist,1]-min(100,stepDelay*int(nrFail[ist]/20))#/stepFraction
                    iD2=prevsol[ist,1]+min(100,stepDelay*(int(nrFail[ist]/20)+1))#/stepFraction
                else:
                    iD1=prevsol[ist,1]
                    iD2=prevsol[ist,1]+stepDelay
                    
                #logging.info("Failure %d : %f %f %f %f %f %f %d "%(ist,iTEC1,iTEC2,stepdTEC,iD1,iD2,stepDelay,nrFail)+str(prevsol[ist]))

             dTECArray=np.arange(iTEC1,iTEC2,stepdTEC)
             dClockArray=np.arange(iD1,iD2,stepDelay)
             datatmp=ph[itm,:,ist]
             #logging.info("getting init par for station %d"%ist)
             if datatmp.count()/float(nF)>0.5:
                 
                 par = getInitPar(datatmp,dTECArray, dClockArray,freq,ClockTECfunc)
                 sol[ist,:]=par[:nparms]
        #wrapflags=np.ones((nSt,nF))     
        for nr_iter in range(2):
            estimate=ClockTECfuncAllStations(freq,sol.T).reshape((nSt,nF)).T
            wraps=np.ma.around(np.divide(estimate-data[itm],2*np.pi))
            data[itm,:]=np.add(2*np.pi*wraps,data[itm])
            #logging.info("fitting masked data itm:%d "%itm + str(data[itm,:].count(axis=0)))
            wrapflags=np.absolute(estimate-data[itm,:])<(1./(nr_iter+1))*np.pi
            #logging.info("flagging dubious wraps" + str(np.sum(np.logical_not(wrapflags),axis=0)))
            #logging.info(str(itm)+":"+str(data[itm,:,-1])+" estimate: "+str(estimate[:,-1])+" "+str(data[itm,:,-1][wrapflags[:,-1]].count())+" "+str(sol[-1]))
            for ist in range(nSt):
                
                if data[itm,:,ist][wrapflags[:,ist]].count()/float(nF)<0.5:
                    logging.info("too many data points flagged t=%d st=%d flags=%d"%(itm,ist,data[itm,:,ist].count()) + str(sol[ist]))
                    sol[ist]=[-10.,-10.]
                    continue
                #print "checking",ist,A.shape,wrapflags.shape,data[itm].shape
                B=A[wrapflags[:,ist]]
                sol[ist]=np.ma.dot(np.linalg.inv(np.dot(B.T,B)),np.ma.dot(B.T,data[itm,:,ist][wrapflags[:,ist]])).T
                #remove jumps in delay
                #if not 'CS' in stations[ist]:
                if initprevsol[ist] and np.abs(np.round((sol[ist,1]-prevsol[ist,1])/steps[1]))>0:
                    sol[ist,:]-=np.round((sol[ist,1]-prevsol[ist,1])/steps[1])*steps
        residual = data[itm] - np.dot(A, sol.T)
        residual = residual - residual[:, 0][:,np.newaxis]
        residual = np.ma.remainder(residual+np.pi, 2*np.pi) - np.pi  
        chi2=np.ma.sum(np.square(np.degrees(residual)),axis=0)/(nF)        

        if returnResiduals:
            residualarray[itm]=residual
            
        chi2select=np.logical_or(np.array(chi2>chi2cut),sol[:,0]<-5)
        if np.any(chi2select):
            logging.info("high chi2 of fit, itm: %d %d "%(itm,np.sum(chi2select)) + str(sol[chi2select])+"stations:" + str(np.arange(nSt)[chi2select])+" chi2 "+str(chi2[chi2select]))
            succes=False
            nrFail[chi2select]+=1
            nrFail[~chi2select]=0
            prevsol[~chi2select]=sol[~chi2select]
            initprevsol[~chi2select]=True
        else:
            prevsol=np.copy(sol)
            succes=True
            initprevsol=np.ones(nSt,dtype=bool)
            nrFail=np.zeros(sol.shape[0],dtype=int)
        tecarray[itm]=sol[:,0]
        clockarray[itm]=sol[:,1]
        if nparms>2:
            fitoffsetarray[itm]=sol[:,2]
    if nparms>2:
        if returnResiduals:
            return tecarray,clockarray,residualarray,fitoffsetarray
        return tecarray,clockarray,fitoffsetarray
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
    flags=tecarray<-5
    resflags=np.logical_or(flags[:,np.newaxis],residualarray==0)
    maskedresiduals=np.ma.array(residualarray,mask=resflags)
    #avgResiduals=np.average(residualarray,axis=0)
    avgResiduals=np.ma.average(maskedresiduals,axis=0)
    wraps,steps=getResidualPhaseWraps(avgResiduals,freq)
    wraps=np.round(wraps-wraps[0])
    logging.info("wraps from residuals: "+str(wraps))
    lats=np.degrees(np.arctan2(pos[:,2],np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1])))
    lats-=lats[0]
    lons=np.degrees(np.arctan2(pos[:,1],pos[:,0]))
    lons-=lons[0]
    lonlat=np.concatenate((lons,lats)).reshape((2,)+lons.shape)
    for nr_iter in range(2):
        TEC=tecarray-tecarray[:,[0]]+steps[0]*(np.round(wraps)-np.round(wraps[0]))
        TEC=np.ma.array(TEC,mask=flags)
        slope=np.ma.dot(np.linalg.inv(np.dot(lonlat,lonlat.T)),np.ma.dot(lonlat,TEC.T))
        chi2=np.ma.sum(np.square(TEC-np.ma.dot(lonlat.T,slope).T),axis=1)/nSt

        chi2select=chi2<np.ma.average(chi2)
        chi2select=chi2<np.ma.average(chi2[chi2select])
        offsets=-1*(np.ma.average(TEC[chi2select]-np.ma.dot(slope.T,lonlat)[chi2select],axis=0))*2.*np.pi/steps[0]
        remainingwraps=np.round(offsets/(2*np.pi))#-np.round(wraps[stationIndices])
        logging.info("offsets: "+str(offsets))
        logging.info("avgTEC: "+str(np.ma.average(TEC[chi2select],axis=0)))
        logging.info("remaining: "+str(remainingwraps))
        wraps+=remainingwraps

        if np.sum(np.absolute(remainingwraps))==0:
            break

    return offsets,wraps,steps     


def doFit(phases,mask,freqs,stations,station_positions,axes,refstIdx='superterp',stationSelect='BA',flagBadChannels=True,flagcut=1.5,chi2cut=30000.,removePhaseWraps=True,combine_pol=False,ignore_stations=["NOTHING_TO_IGNORE",],doFitoffset=False,initSol=[],initoffsets=[]):
    #make sure order of axes is as expected
    stidx=axes.index('ant')
    freqidx=axes.index('freq')
    timeidx=axes.index('time')
    polidx=axes.index('pol')
    data=ma.array(phases,mask=mask).transpose((timeidx,freqidx,stidx,polidx))
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
    refdata=np.ma.sum(np.cos(data[:,:,refstIdx,:])+1j*(np.sin(data[:,:,refstIdx,:])),axis=2)
    refdata=np.arctan2(np.imag(refdata),np.real(refdata))
    #unwrap around mean
    mymean=np.ma.average(refdata,axis=0)
    refdata=np.remainder(refdata-mymean[np.newaxis]+np.pi,2*np.pi)+mymean[np.newaxis]-np.pi
    data-=refdata[:,:,np.newaxis]
    #flag bad channels
    indices=np.arange(nF)
    if flagBadChannels:
        freqselect=np.ones((nF,),dtype=np.bool)
        for nr_iter in range(2):
            rms=np.std(np.std(refdata,axis=0),axis=1)
            freqselect=rms<flagcut*np.average(rms)
            logging.info("iter %d: flagging %d channels"%(nr_iter,np.sum(np.logical_not(freqselect))))
            freqs=freqs[freqselect]
            data=data[:,freqselect]
            refdata=refdata[:,freqselect]
            logging.info("flagging: "+str(indices[np.logical_not(freqselect)]))
            indices=indices[freqselect]
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
    station_positions=station_positions[stationIndices]
    RSstations=[i for i,j in enumerate(stations) if ('RS' in  j)]
    CSstations=[i for i,j in enumerate(stations) if ('CS' in  j)]
    otherstations=[i for i,j in enumerate(stations) if (not ('RS' in j) and not ('CS' in j))]
    logging.info("station indices: "+str(stationIndices)+" RS "+str(RSstations))
    nSt=data.shape[2]
    logging.info("%d selected stations: "%(nSt)+str(stations))   
    # combine polarizationsif requested
    if combine_pol:
        if  npol==2:
            cdata=np.cos(data)+1j*np.sin(data)
            data=np.ma.sum(cdata,axis=3).reshape((nT,nF,nSt,1))
            data=np.arctan2(np.imag(data),np.real(data)) #np.angle doesnot yet return masked array!!
            npol=1
    # guess clock, remove from data
    if len(initSol)<1:
        initclock=getInitClock(data[nT/2:nT/2+100][:,:,RSstations+otherstations],freqs)
        logging.info("initial clocks: "+str(initclock[1])) 
        #init CS clocks to 0

        #logging.info("data before init clock" + str(data[nT/2,:,-1]))
        data[:,:,RSstations+otherstations]=data[:,:,RSstations+otherstations]-(freqs[np.newaxis,:,np.newaxis,np.newaxis]*initclock[1][np.newaxis,np.newaxis,:]*(-1e-9)*2*np.pi)
        #logging.info("clock correction" + str(np.remainder(freqs*initclock[1][-1]*-1e-9*2*np.pi+np.pi,2*np.pi)-np.pi))
        #logging.info("data after init clock" + str(np.remainder(data[nT/2,:,-1]+np.pi,2*np.pi)-np.pi))
    offset=np.zeros((nSt,npol),dtype=np.float32)
    if len(initoffsets)>0:
        offset=initoffsets
        data[:,:,:,:]+=offset[:][np.newaxis,np.newaxis]
    #initialize arrays
    clock=np.zeros((nT,nSt,npol),dtype=np.float32)
    
    tec=np.zeros((nT,nSt,npol),dtype=np.float32)
    if doFitoffset:
        fitoffset=np.zeros((nT,nSt,npol),dtype=np.float32)
    for pol in range(npol):
        #get a good guesss without offset
        #logging.info("sending masked data "+str(data[:,:,:,pol].count()))
        initialchi2cut=chi2cut
        if removePhaseWraps:
            initialchi2cut=1000.
        tecarray,clockarray,residualarray=getClockTECFit(np.ma.copy(data[:,:,:,pol]),freqs,stations,initSol=initSol,returnResiduals=True,chi2cut=initialchi2cut)
        if removePhaseWraps:
            #correctfrist times only,try to make init correct ?
            offset[:,pol],wraps,steps =correctWraps(tecarray[:10000],residualarray[:10000],freqs,station_positions)
        logging.info("residual iter 1, pol %d: "%(pol)+str(residualarray[0,0]))
        logging.info("tec iter 1, pol %d: "%(pol)+str(tecarray[0]))
        logging.info("clock iter 1, pol %d: "%(pol)+str(clockarray[0]))
        if removePhaseWraps:
            logging.info("wraps: "+str(wraps))
        logging.info("offsets: "+str(offset[:,pol]))
        if len(initoffsets)>0:
            offset[:,pol]-=initoffsets[:,pol]
        data[:,:,:,pol]+=offset[:,pol][np.newaxis,np.newaxis]
        if doFitoffset:
            initsol=np.zeros((nSt,3),dtype=np.float32)
            initsol[:,0]=tecarray[0,:]+wraps*steps[0]
            initsol[:,1]=clockarray[0,:]+wraps*steps[1]
            tecarray=0
            clockarray=0
            residualarray=0

            tec[:,:,pol],clock[:,:,pol],fitoffset[:,:,pol]=getClockTECFit(np.ma.copy(data[:,:,:,pol]),freqs,stations,initSol=initsol,returnResiduals=False,chi2cut=chi2cut,nparms=3)
        else:
            if removePhaseWraps:
                initsol=np.zeros((nSt,2),dtype=np.float32)
                initsol[:,0]=tecarray[0,:]+wraps*steps[0]
                initsol[:,1]=clockarray[0,:]+wraps*steps[1]
                logging.info("initsol tec, pol %d: "%(pol)+str(initsol[:,0]))
                logging.info("initsol clock, pol %d: "%(pol)+str(initsol[:,1]))
                tecarray=0
                clockarray=0
                residualarray=0
                tec[:,:,pol],clock[:,:,pol]=getClockTECFit(np.ma.copy(data[:,:,:,pol]),freqs,stations,initSol=initsol,returnResiduals=False,chi2cut=chi2cut)
            else:
                tec[:,:,pol]=tecarray[:,:]
                clock[:,:,pol]=clockarray[:,:]
        logging.info("tec iter 2, pol %d: "%(pol)+str(tec[0,:,pol]))
        logging.info("clock iter 2, pol %d: "%(pol)+str(clock[0,:,pol]))
    if len(initSol)<1:
        clock[:,RSstations+otherstations]+=initclock[1][np.newaxis,:,:]
    if doFitoffset:
        return clock,tec,offset,fitoffset,stations
    return clock,tec,offset,stations
