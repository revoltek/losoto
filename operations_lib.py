#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Some utilities for operations

import sys
import logging
from h5parm import solFetcher


def getParSolsets( step, parset, H ):
    """
    Return the solution-set list for this parset.
    The order is:
    * local step (from the Soltab parameter)
    * if nothing found, global value (from the Soltab parameter)
    * restrict to global Solset var
    * if nothing found use the global Solset var
    * default = all
    """
    allSolsets = H.getSolsets().keys()

    # local val from soltab
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Soltab" ] )
    soltabs = parset.getStringVector( stepOptName, [] )
    solsets = []
    for soltab in soltabs:
        solsets.append(soltab.split('/')[0])

    # global value from soltab
    if solsets == []:
        for soltab in parset.getStringVector( "LoSoTo.Soltab", [] ):
            solsets.append(soltab.split('/')[0])

    # global value from solset
    globalSolsets = parset.getStringVector( "LoSoTo.Solset", allSolsets )
    if solsets != []:
        solsets = list(set(solsets).intersection(set(globalSolsets)))
    else:
        solsets = globalSolsets

    # default value
    if solsets == []:
        solsets = allSolsets

    # sanity check
    for solset in solsets[:]:
        if solset not in allSolsets:
            logging.warning("Solution-set \""+solset+"\" not found. Ignoring")
            solsets.remove(solset)

    return list(set(solsets))


def getParSoltabs( step, parset, H ):
    """
    Return the solution-table list (in ["solset/soltab",...] form) for this step.
        - compatible with the solset parameters
        - compatible the soltype parameters
    The order is:
    * local step value
    * global value
    * default = all
    """

    # local val
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Soltab" ] )
    ssst = parset.getStringVector( stepOptName, [] )

    # global value
    if ssst == []:
        ssst = parset.getStringVector( "LoSoTo.Soltab", [] )

    # default value
    if ssst == []:
        # add all the table in the available Solsets
        for solset in getParSolsets( step, parset, H ):
            for soltab in H.getSoltabs(solset).keys():
                ssst.append(solset+'/'+soltab)

    # sanity check
    allawedSolTypes = getParSolTypes( step, parset, H )
    for s in ssst[:]:
        solset, soltab = s.split('/')
        # check that soltab exists and that the declared solset is usable
        if solset not in getParSolsets( step, parset, H ) or \
                soltab not in H.getSoltabs(solset).keys():
            logging.warning("Solution-table \""+ solset+"/"+soltab+"\" not found. Ignoring.")
            ssst.remove(s)
        # check if the soltab is compatible with the chosen solTypes
        elif H.getSoltab(solset, soltab)._v_title not in allawedSolTypes and allawedSolTypes != []:
            ssst.remove(s)

    return ssst


def getParSolTypes( step, parset, H ):
    """
    Return the SolType list for this step.
    The order is:
    * local step value
    * global value
    * default = [] (==all)
    """

    # local val
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "SolType" ] )
    solTypes = parset.getStringVector( stepOptName, [] )

    # global val or default
    if solTypes == []:
        solTypes = parset.getStringVector( "LoSoTo.SolType", [] )

    return solTypes


def getParAxis( step, parset, H, axisName ):
    """
    Return the axis val array for this step.
        - check if all the soltabs have this axis.
    The order is:
    * local step value
    * global value
    * default = None (which keep all in setSelection)
    """
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, axisName.lower() ] )
    # local val
    axisVals = parset.getString( stepOptName, '' )
    # if the user defined a vector, load it as a vector, otherwise keep string
    if axisVals != '' and xisVals[0] == '[' and axisVals[-1] == ']':
        axisVals = parset.getStringVector( stepOptName, [] )

    # global val
    if axisVals == [] or axisVals == '':
        axisVals = parset.getString( "LoSoTo."+axisName.lower(), '' )
        # if the user defined a vector, load it as a vector, otherwise keep string
        if axisVals != '' and axisVals[0] == '[' and axisVals[-1] == ']':
            axisVals = parset.getStringVector( "LoSoTo."+axisName.lower(), [] )

    # default val
    if axisVals == [] or axisVals == '': axisVals = None

    return axisVals


def openSoltabs( H, ss_sts ):
    """
    Return a list of soltab objects
    Keyword arguments:
    ss_sts -- 'solution-set/solution-tabs' list
    """
    soltabs = []
    for ss_st in ss_sts:
        ss, st = ss_st.split('/')
        soltabs.append( H.getSoltab(ss, st) )

    return soltabs

def removeKeys( dic, keys = [] ):
    """
    Remove a list of keys from a dict and return a new one.
    Keyword arguments:
    dic -- the input dictionary
    keys -- a list of arguments to remove or a string for single argument
    """
    dicCopy = dict(dic)
    if type(keys) is str: keys = [keys]
    for key in keys:
        del dicCopy[key]
    return dicCopy


def unwrap_fft(phase, iterations=3):
    """
    Unwrap phase using Fourier techniques.

    For details, see:
    Marvin A. Schofield & Yimei Zhu, Optics Letters, 28, 14 (2003)

    Keyword arguments:
    phase -- array of phase solutions
    iterations -- number of iterations to perform
    """
    import numpy as np

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


def unwrap_huib( x, window = 10, alpha = 0.01, iterations = 3,
    clip_range = [ 170., 180. ] ):
    """
    Unwrap the x array, if it is shorter than 2*window, use np.unwrap()
    """

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
