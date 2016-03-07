import matplotlib 
matplotlib.use('agg') 
import numpy as np
import pylab
import time
import multiprocessing
import os
import pyrap.quanta as qa
import pyrap.measures as pm
import ephem
import rad2hmsdms

import numpy as np
import time
import pylab

from reformat import *
import getpass
import sys

import tables
from progressbar import ProgressBar
import scipy.ndimage

class ClassMovie():
    def __init__(self,Hdf5File,TStart=0,TEnd=-1,DoRM=True,incr=1):

        self.Hdf5File=Hdf5File

        H=tables.open_file(self.Hdf5File)



        self.pngDir="png"
        self.incr=incr

        os.system("rm -rf %s"%self.pngDir)
        os.system("mkdir -p %s"%self.pngDir)
        #if not(os.path.exists(self.pngDir)): os.makedirs(self.pngDir)
        
        self.AntNames=H.root.sol000.antenna[:]["name"]
        self.Freqs=H.root.sol000.amplitude000.freq[:]/1.e6
        self.Times=H.root.sol000.amplitude000.time[:]

        self.amplitudes=H.root.sol000.amplitude000.val[:]
        self.phase=H.root.sol000.phase000.val[:]


        self.na=self.AntNames.size
        self.NTimes=self.Times.size
        self.NFreqs=self.Freqs.size
        ra,dec=H.root.sol000.source[:]["dir"][0]
        self.decdir=dec
        self.radir=ra
        self.TStart=TStart
        self.TEnd=TEnd
        if TEnd==-1:
            self.TEnd=self.NTimes-1
        self.DoOff=False

    def test(self,t=0):
        J=self.GiveJones(t)
        fig=pylab.figure(1)
        self.setCoord()
        self.plot_individual2(J,t,fig)
        pylab.draw()
        pylab.show()

    def GiveJones(self,t):
        ind=np.argmin(np.abs(self.Times-t))

        amp=self.amplitudes[:,0,:,:,ind]

        pha=self.phase[:,0,:,:,ind]
        Jones=amp*np.exp(2.*np.pi*1j*pha)
        # x/y, ant, freq
        Jones=np.swapaxes(Jones,0,1)
        Jones=np.swapaxes(Jones,1,2)
        Jout=np.zeros((self.na,self.NFreqs,2,2),dtype=Jones.dtype)
        Jout[:,:,0,0]=Jones[:,:,0]
        Jout[:,:,1,1]=Jones[:,:,1]
        

        return Jout

    def setCoord(self):
        na=self.na
        ny=int(np.sqrt(na/1.77))
        nx=na/ny
        ny+=1

        xcoord=np.linspace(0.05,.95,nx+1)
        ycoord=np.linspace(0.05,.95,ny+1)
    
        self.acoord=np.zeros((na,2),dtype=np.float)
        for ant in range(na)[::-1]:
            i=int(ant/nx)
            j=ant-i*nx
            self.acoord[ant,0]=xcoord[j]
            self.acoord[ant,1]=ycoord[i]

        self.dx=xcoord[1]-xcoord[0]
        self.dy=(ycoord[1]-ycoord[0])/2.
        self.margx=4e-3
        self.margy=2.e-2

        ArrMed=np.median(np.median(np.abs(self.amplitudes[0]),axis=1),axis=0)
        ArrMed=np.max(scipy.ndimage.median_filter(ArrMed,5))
        
        self.ylim=[0.,np.max(ArrMed)*1.2]
        self.xlim=[np.min(self.Freqs),np.max(self.Freqs)]

        
    def plot_parallel(self):
        
        self.setCoord()
        AntNames=self.AntNames
        Freqs=self.Freqs
        Times=self.Times

        na=self.na
        NTimes=self.NTimes
        NFreqs=self.NFreqs
        decdir=self.decdir
        radir=self.radir

        NCPU=20
        self.figs=[]
        for i in range(NCPU):
            self.figs.append(pylab.figure(i,figsize=(12,7)))
        #global acoord, margx, margy, dx, dy, xlim, ylim
    
        t0=Times[0]
    

        jobs=[]
        ii=0
        pylab.clf()
        #work_queue = multiprocessing.JoinableQueue()
        work_queue = multiprocessing.Queue()
        for t in range(self.TStart,self.TEnd,self.incr):#NTimes):
            J=self.GiveJones(self.Times[t])
            jobs.append([J.copy(),t])
            #work_queue.put([J.copy(),t])


        work_queue = multiprocessing.JoinableQueue()
        result_queue = multiprocessing.JoinableQueue()

        for job in jobs:
            work_queue.put(job)



        workerlist=[]
        for ii in range(NCPU):
            workerlist.append(WorkerPlotGains(work_queue, result_queue,self))
            workerlist[ii].start()
    
        results = []
        lold=len(results)

        ss=jobs
        lold=len(results)
        
        toolbar_width = 50
        pBAR= ProgressBar('white', block='=', empty=' ',Title="Solving")
        pBAR.render(0, '%i/%i' % (0,len(ss)-1.))

        while len(results) < len(jobs):
            result = result_queue.get()
            results.append(result)

            #print result
            if len(results)>lold:
                lold=len(results)
                pBAR.render(int(100* float(lold) / (len(ss)-1.)), 'step %i/%i' % (lold,len(ss)-1.))

        # del(workerlist)
        # del(results)


        # import errno
        # for ii in range(NCPU):
        #     notintr = False
        #     while not notintr:
        #         try:
        #             workerlist[i].join() # "Offending code"
        #             notintr = True
        #         except OSError, ose:
        #             if ose.errno != errno.EINTR:
        #                 raise ose

        for ii in range(NCPU):
            workerlist[ii].shutdown()
            workerlist[ii].terminate()
            workerlist[ii].join()
        del(workerlist)

        # sys.exit()


    def plot_individual2(self,Jones,t,fig):

        AntNames=self.AntNames
        Freqs=self.Freqs
        Times=self.Times

        na=self.na
        NTimes=self.NTimes
        NFreqs=self.NFreqs
        decdir=self.decdir
        radir=self.radir
        acoord=self.acoord
        dx=self.dx
        dy=self.dy
        margx=self.margx
        margy=self.margy
        ylim=self.ylim
        xlim=self.xlim


        tt=Times[t]-Times[0]
        th=int(tt/3600.)
        tt=tt-th*3600.
        tm=int(tt/60.)
        ts=tt-tm*60.
        #srttime="Since beginning: %2.2ih%2.2im%05.2fs"%(th,tm,ts)
        for ant in range(na):
            a = fig.add_axes([acoord[ant,0]+margx, acoord[ant,1]+margy, dx-2.*margx,dy-margy])#, axisbg='y')
            if self.DoOff:
                a.plot(Freqs,np.abs(Jones[ant,:,0,1]),marker='.',ls='',mew=0.1,ms=1.5,color="gray")
                a.plot(Freqs,np.abs(Jones[ant,:,1,0]),marker='.',ls='',mew=0.1,ms=1.5,color="gray")
            a.plot(Freqs,np.abs(Jones[ant,:,0,0]),marker='.',ls='',mew=0.1,ms=1.5,color="black")
            a.plot(Freqs,np.abs(Jones[ant,:,1,1]),marker='.',ls='',mew=0.1,ms=1.5,color="black")
            if ant==0:
                a.set_xlabel("Freq")
                a.set_ylabel("Amp")
            pylab.setp(a, xticks=[], yticks=[],ylim=self.ylim)#,xlim=xlim
            a = fig.add_axes([acoord[ant,0]+margx, acoord[ant,1]+dy, dx-2.*margx,dy-margy])#, axisbg='y')
            if self.DoOff:
                a.plot(Freqs,np.angle(Jones[ant,:,0,1])-np.angle(Jones[2,:,0,1]),marker='.',ls='',mew=0.1,ms=1.5,color="gray")
                a.plot(Freqs,np.angle(Jones[ant,:,1,0])-np.angle(Jones[2,:,1,0]),marker='.',ls='',mew=0.1,ms=1.5,color="gray")
            a.plot(Freqs,np.angle(Jones[ant,:,0,0])-np.angle(Jones[2,:,0,0]),marker='.',ls='',mew=0.1,ms=1.5,color="black")
            a.plot(Freqs,np.angle(Jones[ant,:,1,1])-np.angle(Jones[2,:,1,1]),marker='.',ls='',mew=0.1,ms=1.5,color="black")
            a.title.set_text(AntNames[ant])
            a.title.set_fontsize(8)
            pylab.setp(a, xticks=[], yticks=[],ylim=[-np.pi,np.pi])
            
            if ant==0:
                a.set_ylabel("Pha")
            #a.annotate(str(ant), xy=(acoord[ant,0]+(dx-2.*marg)/2., acoord[ant,1]+marg),  xycoords='figure fraction')
            sumtrace=0.
            sumdet=0.
            for ch in range(NFreqs):
                sumtrace+=np.trace(Jones[ant,ch,:,:])
                sumdet+=np.linalg.det(Jones[ant,ch,:,:])
        sumtrace/=NFreqs
        sumdet/=NFreqs
        #a.annotate("%s (%4.1f,%4.1f)"%(AntNames[ant],sumtrace,sumdet), xy=(acoord[ant,0], acoord[ant,1]+dy+dy-margy),  xycoords='figure fraction',fontsize=8)
        #a.annotate("%s\n%6.1f,%6.1f"%(AntNames[ant],np.log10(np.abs(sumtrace)),np.log10(np.abs(sumdet))), xy=(acoord[ant,0]+dx/2.,\
        #                            acoord[ant,1]+dy+dy-margy),  xycoords='figure fraction',fontsize=8, horizontalalignment="center")
    
        ttt=Times[t]
        strd,az,alt=self.give_elevation(ttt)
        srttime="%s, (az, alt)=(%05.2f, %05.2f) deg, Freqs=[%5.1f, %5.1f]MHz, log Tr. and log det. below Name"%(strd,az*180./np.pi,alt*180./np.pi,np.min(Freqs),np.max(Freqs))
        
        a.annotate(srttime, xy=(0.5, 0.97),  xycoords='figure fraction', horizontalalignment="center")
    
    def give_elevation(self,tt):
    
        AntNames=self.AntNames
        Freqs=self.Freqs
        Times=self.Times

        na=self.na
        NTimes=self.NTimes
        NFreqs=self.NFreqs
        decdir=self.decdir
        radir=self.radir

        lofar = ephem.Observer()
        lofar.lat, lofar.long = '52.834444', '6.370556'
        time_center = qa.quantity(tt, 's')
        me = pm.measures()
        dict_time_center_MDJ = me.epoch('utc', time_center)
        time_center_MDJ=dict_time_center_MDJ['m0']['value']
        JD=time_center_MDJ+2400000.5-2415020
        d=ephem.Date(JD)
        lofar.date=d
        y,mth,d,h,m,s=d.tuple()
        strdate="%2.2i/%2.2i/%4.4i %2.2i:%2.2i:%2.2i"%(d,mth,y,h,m,s)
        rastr =rad2hmsdms.rad2hmsdms(radir,Type="ra").replace(" ",":")
        decstr=rad2hmsdms.rad2hmsdms(decdir,Type="dec").replace(" ",":")
    
        db="Bootes,f|M|F7,"+rastr+","+decstr+",2.02,2000"
        o=ephem.readdb(db)
        o.compute(lofar)
        az,alt= o.az.__float__(),o.alt.__float__()
    
        return strdate,az,alt
    
    
                
    
    
    
    def make_movie(self,name):
        ss="mencoder -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:vbitrate=6160000:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:"+\
            "mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq -mf type=png:fps=20 -nosound -o "+name+".mpg mf://\png/*.png  > lala_mpg 2>&1"
        os.system(ss)
    
class WorkerPlotGains(multiprocessing.Process):
    def __init__(self,
            work_queue,
            result_queue,ClassMovie):
        multiprocessing.Process.__init__(self)
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.exit = multiprocessing.Event()
        self.CMovie=ClassMovie

    def shutdown(self):
        self.exit.set()
    def run(self):
        while not self.kill_received:
            try:
                job = self.work_queue.get()
            except:
                break

            Jin,t=job

            ii=multiprocessing.current_process()._identity[0]-1
            self.CMovie.figs[ii].clf()
            self.CMovie.plot_individual2(Jin,t,self.CMovie.figs[ii])
            self.CMovie.figs[ii].savefig("png/lala%4.4i.png"%t)


            self.result_queue.put([1])



if __name__=="__main__":
    import sys
    Obs=sys.argv[1]
    if len(sys.argv)==2:
        TStart=0
        TEnd=-1
    else:
        TStart=int(sys.argv[2])
        TEnd=int(sys.argv[3])

    incr=int(sys.argv[4])
    PickleFile=sys.argv[5]
    import ModPickle
    Setup=ModPickle.Load(PickleFile)

    M=ClassMovie(Setup,TStart=TStart,TEnd=TEnd,DoRM=False,incr=incr)
    M.plot_parallel()

        
    #class ClassMovie():
    #    def __init__(self,MainDir="/home/%s/PipeSurvey/",Obs="L74464",TStart=0,TEnd=-1):
