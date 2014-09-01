#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Some utilities for operations

import sys
import logging
from h5parm import solFetcher

#def getParAnts( step, parset, H ):
#    """
#    Return the Ant array for this step.
#    If more then one solset is involved, return the intersection
#    The order is:
#    * local step value
#    * global value
#    * default = all
#    """
#    allAnts = []
#    for solset in getParSolsets(step, parset, H):
#        allAnts.append(set(H.getAnt(solset).keys()))
#    allAnts = list(set.intersection(*allAnts))
#
#    # local val
#    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Ant" ] )
#    ants = parset.getStringVector( stepOptName, [] )
#
#    # global val or default
#    if ants == []:
#        ants = parset.getStringVector( "LoSoTo.Ant", allAnts )
#
#    # sanity check
#    for ant in ants[:]:
#        if ant not in allAnts:
#            logging.warning("Antenna \""+ant+"\" not found. Ignoring.")
#            ants.remove(ant)
#
#    return ants


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


#def getParDirs( step, parset, H ):
#    """
#    Return the directions array for this step.
#        - check if at least one solset has this direction.
#    The order is:
#    * local step value
#    * global value
#    * default = []
#    """
#    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Dir" ] )
#    # local val
#    dirs = parset.getStringVector( stepOptName, [] )
#
#    # global val or default
#    if dirs == []:
#        dirs = parset.getStringVector( "LoSoTo.Dir", [] )
#
#    # check that directions are valid for at least one solset in step
#    solsets = getParSolsets( step, parset, H )
#    for dir in dirs[:]:
#        found = False
#        for solset in solsets:
#            if dir in H.getSou(solset).keys():
#                found = True
#                break
#        if not found:
#            logging.warning("Direction \""+dir+"\" not found. Ignoring.")
#            dirs.remove(dir)
#
#    return dirs


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
    axisVals = parset.getStringVector( stepOptName, [] )

    # global val
    if axisVals == []:
        axisVals = parset.getStringVector( "LoSoTo."+axisName.lower(), [] )

    # default val
    if axisVals == []: axisVals = None

    # check that pols are valid for at least one soltab in step
#    solsets = getParSolsets( step, parset, H )
#    for axisVal in axisVals[:]:
#        found = False
#        for solset in solsets:
#            soltabs = H.getSoltabs(solset=solset)
#            for soltab_name in soltabs.keys():
#                sf = solFetcher(soltabs[soltab_name])
#                if axisVal in sf.getValuesAxis(axisName):
#                    found = True
#                    break
#            if found:
#                break
#        if not found:
#            logging.warning("Values \""+axisVal+"\" not found. Ignoring.")
#            axisVals.remove(axisVal)

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


def unwrap_fft(phase, iterations=1):
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
