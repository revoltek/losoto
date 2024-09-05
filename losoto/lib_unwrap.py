#!/usr/bin/python
# -*- coding: utf-8 -*-

#from __future__ import division

import itertools, math
import numpy as np
import scipy.fftpack as fft
from losoto._logging import logger as logging

def unwrap_fft(phase, iterations=3):
    """
    Unwrap phase using Fourier techniques.

    For details, see:
    Marvin A. Schofield & Yimei Zhu, Optics Letters, 28, 14 (2003)

    Keyword arguments:
    phase -- array of phase solutions
    iterations -- number of iterations to perform
    """
    puRadius=lambda x : np.roll( np.roll(
          np.add.outer( np.arange(-x.shape[0]/2+1,x.shape[0]/2+1)**2.0,
                        np.arange(-x.shape[1]/2+1,x.shape[1]/2+1)**2.0 ),
          x.shape[1]/2+1,axis=1), x.shape[0]/2+1,axis=0)+1e-9

    idt,dt=np.fft.ifft2,np.fft.fft2
    puOp=lambda x : idt( np.where(puRadius(x)==1e-9,1,puRadius(x)**-1.0)*dt(
          np.cos(x)*idt(puRadius(x)*dt(np.sin(x)))
         -np.sin(x)*idt(puRadius(x)*dt(np.cos(x))) ) )

    def phaseUnwrapper(ip):
       mirrored=np.zeros([x*2 for x in ip.shape])
       mirrored[:ip.shape[0],:ip.shape[1]]=ip
       mirrored[ip.shape[0]:,:ip.shape[1]]=ip[::-1,:]
       mirrored[ip.shape[0]:,ip.shape[1]:]=ip[::-1,::-1]
       mirrored[:ip.shape[0],ip.shape[1]:]=ip[:,::-1]

       return (ip+2*np.pi*
             np.round((puOp(mirrored).real[:ip.shape[0],:ip.shape[1]]-ip)
             /2/np.pi))

    phase2D = phase[:, None]
    i = 0
    if iterations < 1:
        interations = 1
    while i < iterations:
        i += 1
        phase2D = phaseUnwrapper(phase2D)

    return phase2D[:, 0]


def unwrap(phase, window_size=5):
    """
    Unwrap phase by estimating the trend of the phase signal.
    """
    # Allocate result.
    out = np.zeros(phase.shape)

    windowl = np.array([math.fmod(phase[0], 2.0 * math.pi)] * window_size)

    delta = math.fmod(phase[1] - windowl[0], 2.0 * math.pi)
    if delta < -math.pi:
        delta += 2.0 * math.pi
    elif delta > math.pi:
        delta -= 2.0 * math.pi
    windowu = np.array([windowl[0] + delta] * window_size)

    out[0] = windowl[0]
    out[1] = windowu[0]

    meanl = windowl.mean()
    meanu = windowu.mean()
    slope = (meanu - meanl) / float(window_size)

    for i in range(2, len(phase)):
        ref = meanu + (1.0 + (float(window_size) - 1.0) / 2.0) * slope
        delta = math.fmod(phase[i] - ref, 2.0 * math.pi)

        if delta < -math.pi:
            delta += 2.0 * math.pi
        elif delta > math.pi:
            delta -= 2.0 * math.pi

        out[i] = ref + delta

        windowl[:-1] = windowl[1:]
        windowl[-1] = windowu[0]
        windowu[:-1] = windowu[1:]
        windowu[-1] = out[i]

        meanl = windowl.mean()
        meanu = windowu.mean()
        slope = (meanu - meanl) / float(window_size)

    return out


def unwrap_huib( x, window = 10, alpha = 0.01, iterations = 3,
    clip_range = [ 170., 180. ] ):
    """
    Unwrap the x array, if it is shorter than 2*window, use np.unwrap()
    """
    import np as np

    if len(x) < 2*window: return np.unwrap(x)

    xx = np.array( x, dtype = np.float64 )
#   a = zeros( ( window ), dtype = np.float64 )
#   a[ -1 ] = 1.
    if ( len( clip_range ) == 2 ):
        o = clip_range[ 0 ]
        s = ( clip_range[ 1 ] - clip_range[ 0 ] ) / 90.
    a = np.ones( ( window ), dtype = np.float64 ) / float( window )
    xs = xx[ 0 ]
    for j in range( 2 * iterations ):
        for k in range( window, len( x ) ):
            xi = xx[ k - window : k ]
            xp = np.dot( xi, a )
            e = xx[ k ] - xp
            e = np.mod( e + 180., 360. ) - 180.
            if ( len( clip_range ) == 2 ):
                if ( abs( e ) > o ):
                    e = sign( e ) * ( s * degrees( atan( radians( ( abs( e ) - o ) / s ) ) ) + o )
            xx[ k ] = xp + e
#           a = a + xx[ k - window : k ] * alpha * e / 360.
            a = a + xi * alpha * e / ( np.dot( xi, xi ) + 1.e-6 )
        xx = xx[ : : -1 ].copy()
    xx = xx - xx[ 0 ] + xs
    return xx

# This script implements a DCT-based one-step path independent phase
# unwrapping method described in DOI: 10.3390/jimaging1010031 .
# Phase will be unwrapped up to some additive constant.
# - Unwrapping performance depends on the phase map's noise level and C2 smoothness.
# - If zeropadding has to be applied beforehand for FFTs to work, some ugly
#   boundary effects might occur. Best use dims supported by FFTs.
# - Masking out parts of the phase map with zeros also leads to boundary effects.
# - Use it stand alone or as initial guess for preconditioned conjugate gradient method.

# 2D DCT and inverse
def dct2(arr, inverse=False):
    if inverse:
        return np.transpose(fft.idct(np.transpose(fft.idct(arr,norm='ortho')),norm='ortho'))
    else:
        return np.transpose(fft.dct(np.transpose(fft.dct(arr,norm='ortho')),norm='ortho'))

# Laplace operator and inverse
def laplacian(arr, inverse=False):
    # precompute Fourier domain coordinates for the Laplacians of phase_wrapped
    m = np.arange(arr.shape[1])
    n = np.arange(arr.shape[0])
    m,n = np.meshgrid(m,n)
    m = m.astype(float)
    n = n.astype(float)
    m[:,0] = 1/np.sqrt(arr.shape[1])
    n[0,:] = 1/np.sqrt(arr.shape[0])
    coord = m*m+n*n
    m = n = None

    if inverse:
        return dct2(dct2(arr)/coord, inverse=True) # factor -M*N/4/pi**2 omitted
    else:
        return dct2(dct2(arr)*coord, inverse=True) # factor -4*pi**2/M/N omitted

# phase unwrapping, see DOI: 10.3390/jimaging1010031
def unwrap_2d(arr, flags = None, coord_x = None, coord_y = None):
    """
    if flags are specified do interp
    """
    if not flags is None:
        import scipy.interpolate
        if coord_x is None or coord_y is None:
            logging.error('Cannot unwrap with flags and no coordinates.')
            return
        shapeOrig = arr.shape
        grid = np.array([x for x in itertools.product(coord_x,coord_y)])
        arr = np.ndarray.flatten(arr)
        flags = np.ndarray.flatten(flags)
        arr[flags] = scipy.interpolate.griddata(grid[~flags], arr[~flags], grid[flags], 'nearest')
        arr = arr.reshape(shapeOrig)
    return arr + np.round( ( laplacian( np.cos(arr)*laplacian(np.sin(arr)) - np.sin(arr)*laplacian(np.cos(arr)), inverse=True ) - arr ) / 2/np.pi ) * 2*np.pi

#if __name__ == "__main__":
#
#    import matplotlib.pyplot as pl
#    from time import time
#
#    # phase map dims, zeropadding and masking should be avoided if possible
#    N = 1024
#    M = 1024
#    
#    # make some coordinates for the synthetic phase map
#    x = np.arange(N)
#    y = np.arange(M)
#    x,y = np.meshgrid(x,y)
#    x = x / (np.max(x)/(2*np.pi))
#    y = y / (np.max(y)/(2*np.pi))
#    
#    snr = 1 # SNR for the Gaussian noise
#    
#    # some synthetic phase maps
#    #phase_orig = np.sin(x)*x+np.cos(y)*y
#    phase_orig = np.sin(x*y)*x+np.cos(x*x+y*y)*y
#    #phase_orig = (x-np.max(x)/2)**2 + (y-np.max(y)/2)**2
#    x = y = None
#    phase_orig /= np.max(np.abs(phase_orig)) / (np.pi*4.) # scale phase map amplitude
#    phase_orig += np.array( np.random.normal(scale=np.sqrt(0.5/snr),size=(N,M)).astype(np.float32) ) # add noise
#    
#    phase_wrapped = np.arctan2( np.sin(phase_orig), np.cos(phase_orig) ) # wrap phase map to [-pi,+pi]
#    
#    # precompute Fourier domain coordinates for the Laplacians of phase_wrapped
#    m = np.arange(N)
#    n = np.arange(M)
#    m,n = np.meshgrid(m,n)
#    m = m.astype(float)
#    n = n.astype(float)
#    m[:,0] = (1./np.sqrt(M))
#    n[0,:] = (1./np.sqrt(N))
#    coord = m*m+n*n
#    m = n = None
#    
#    #k = 1 # timing iterations
#    t0=time()
#    phase_unwrapped = unwrap_dct(phase_wrapped) # warm-up
#    #for i in range(k):
#    #    phase_unwrapped = unwrap_dct(phase_wrapped)
#    t1 = time()-t0
#    coord = None
#    #print('time:\t%fs'%(t1))
#    
#    #print('norm:', np.abs(phase_unwrapped-phase_orig) )
#    
#    fig, ((ax1,ax2),(ax3,ax4)) = pl.subplots(2,2)
#    im1 = ax1.imshow(phase_orig.__array__())
#    pl.colorbar(im1,ax=ax1)
#    ax1.set_title('original phase, SNR=%i'%snr)
#    im2 = ax2.imshow(phase_wrapped.__array__())
#    pl.colorbar(im2,ax=ax2)
#    ax2.set_title('wrapped phase')
#    im3 = ax3.imshow(phase_unwrapped.__array__())
#    pl.colorbar(im3,ax=ax3)
#    ax3.set_title('unwrapped phase')
#    im4 = ax4.imshow((phase_unwrapped-phase_orig).__array__())
#    pl.colorbar(im4,ax=ax4)
#    ax4.set_title('difference')
#    pl.savefig('test.png')
