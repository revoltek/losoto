from __future__ import division
import numpy as np
cimport numpy as np
import time

    
cimport cython
cdef extern from "complex.h":
    complex conj(complex z)

@cython.wraparound(False)
@cython.boundscheck(False)
def FastInv22(np.ndarray[complex, ndim=2] A, int H):
    cdef complex ff,a,b,c,d
    if H==0:
        a=A[0, 0]
        b=A[0, 1]
        c=A[1, 0]
        d=A[1, 1]
    else:
        a=conj(A[0, 0])
        b=conj(A[1, 0])
        c=conj(A[0, 1])
        d=conj(A[1, 1])
    ff=1./((a*d-c*b))
    A[0, 0]=ff*d
    A[0, 1]=-ff*b
    A[1, 0]=-ff*c
    A[1, 1]=ff*a


@cython.wraparound(False)
@cython.boundscheck(False)
def Prod22(np.ndarray[complex, ndim=2] A, np.ndarray[complex, ndim=2] B, np.ndarray[complex, ndim=2] Out):
    Out[:,:]=0
    #cdef np.ndarray[complex, ndim=2] Out=np.zeros((2,2), dtype=complex)
    Out[0, 0]=A[0, 0]*B[0, 0]+A[0, 1]*B[1, 0]
    Out[0, 1]=A[0, 0]*B[0, 1]+A[0, 1]*B[1, 1]
    Out[1, 0]=A[1, 0]*B[0, 0]+A[1, 1]*B[1, 0]
    Out[1, 1]=A[1, 0]*B[0, 1]+A[1, 1]*B[1, 1]
    return Out

from progressbar import ProgressBar


#Sols=(na,nt,2,2)
@cython.wraparound(False)
@cython.boundscheck(False)
def FastApply(np.ndarray[complex, ndim=3] Vis, np.ndarray[double, ndim=1] TIME, np.ndarray[int, ndim=1] A0, \
                  np.ndarray[int, ndim=1] A1, np.ndarray[complex, ndim=4] Gains, double t0sols, double DtSol):

    cdef int Nrows = Vis.shape[0]
    cdef int Nchan = Vis.shape[1]
    cdef int NTimeSols = Gains.shape[1]
    cdef unsigned int irow = 0
    cdef unsigned int ichan = 0

    cdef np.ndarray[complex, ndim=2] VisBuf
    cdef np.ndarray[complex, ndim=2] tmp0=np.zeros((2,2), dtype=complex)
    cdef np.ndarray[complex, ndim=2] tmp1=np.zeros((2,2), dtype=complex)
    cdef np.ndarray[complex, ndim=2] tmpSol0=np.zeros((2,2), dtype=complex)
    cdef np.ndarray[complex, ndim=2] tmpSol1=np.zeros((2,2), dtype=complex)

    cdef double time
    cdef unsigned int itime_sol,iA0, iA1
    pBAR= ProgressBar('white', block='=', empty=' ',Title="Apply Jones")

    cdef int icount=0

    for irow from 0 <= irow < Nrows:
        time=TIME[irow]
        itime_sol=int((time-t0sols)/DtSol)
        if itime_sol<0: itime_sol=0
        if itime_sol>(NTimeSols-1): itime_sol=NTimeSols-1
        iA0=A0[irow]
        iA1=A1[irow]
        
        # This:
        # [station, time, freqs, J0, J1]
        #tmpSol0[:,:]=Gains[iA0,itime_sol,:,:]
        #tmpSol1[:,:]=Gains[iA1,itime_sol,:,:]
        # Globalbd:
        # In [6]: H.amplitudes[:].shape
        # Out[6]: (7190, 242, 53, 1, 2, 2)
        tmpSol0[:,:]=Gains[itime_sol,iA0,:,:]
        tmpSol1[:,:]=Gains[itime_sol,iA1,:,:]
        FastInv22(tmpSol0,H=0)
        FastInv22(tmpSol1,H=1)


        for ichan from 0 <= ichan < Nchan:
            VisBuf=Vis[irow,ichan,:].reshape(2,2)
            Prod22(tmpSol0,VisBuf,tmp0)
            Prod22(tmp0,tmpSol1,VisBuf)

        icount+=1

        if (icount==10127)|(irow==Nrows-1):
            pBAR.render(int(100*float(irow+1)/Nrows), '%i/%i' % (irow+1,Nrows))
            icount=0
    
    return Vis
            

