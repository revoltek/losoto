#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Some utilities for operations.py

import sys
import logging

# TODO: add pattern matching
def getParAnts( step, parset, H ):
    """
    Return the Ant array for this step.
    If more then one solset is involved, return the intersection
    The order is:
    * local step value
    * global value
    * default = all
    """
    allAnts = []
    for solset in getParSolsets(step, parset, H):
        allAnts.append(set(H.getAnt(solset).keys()))
    allAnts = list(set.intersection(*allAnts))

    # local val
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Ant" ] )
    ants = parset.getStringVector( stepOptName, [] )
    
    # global val or default
    if ants == []:
        ants = parset.getStringVector( "LoSoTo.Ant", allAnts )

    # sanity check
    for ant in ants:
        if ant not in allAnts:
            logging.error("Cannot find Ant", Ant, ", ignoring.")
            ants.remove(Ant)

    return ants


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
        # add all the table in the globalSolset if defined,
        # otherwise add all tables
        for solset in getParSolsets( step, parset, H ):
            for soltab in H.getSoltabs(solset).keys():
                ssst.append(solset+'/'+soltab)

    # sanity check
    for s in ssst:
        solset, soltab = s.split('/')
        # check that soltab exists and that the declared solset is usable
        if soltab not in H.getSoltabs(solset).keys() or\
                solset not in getParSolsets( step, parset, H ):
            logging.error("Solution-table", soltab, " not available, ignoring.")
            ssst.remove(s)

    return ssst


def getParSolsets( step, parset, H ):
    """
    Return the solution-set list for this parset.
    The order is:
    * local step (from the Soltab parameter)
    * global value (from the Soltab + Solset parameter)
    * default = all
    """
    allSolsets = H.getSolsets().keys()

    # local val
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Soltab" ] )
    soltabs = parset.getStringVector( stepOptName, [] )
    solsets = []
    for soltab in soltabs:
        solsets.append(soltab.split('/')[0])

    # global value from soltab
    if solsets == []:
        for soltab in parset.getStringVector( "LoSoTo.Soltab", [] ):
            solsets.extend(soltab.split('/')[0])

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
    for solset in solsets:
        if solset not in allSolsets:
            logging.error("Solution-set", solset, " in not in the HDF5 file, ignoring")
            solsets.remove(solset)

    return list(set(solsets))


def getParSolType( step, parset, H ):
    """
    Return the SolType list for this step.
    The order is:
    * local step value
    * global value
    * default = [] (==all)
    """

    # local val
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "SolType" ] )
    SolType = parset.getStringVector( stepOptName, [] )
    
    # global val or default
    if SolTypes == []:
        SolTypes = parset.getStringVector( "LoSoTo.solType", [] )

    return solTypes



def getParDirs( step, parset, H ):
    """
    Return the directions array for this step.
        - check is all the toltb has this direction.
    The order is:
    * local step value
    * global value
    * default = []
    """
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Dir" ] )
    # local val
    dirs = parset.getStringVector( stepOptName, [] )
    if dirs == []:
        dirs = parset.getStringVector( "LoSoTo.Dir", [] )
        # global val or default

    return dirs


def getParPols( step, parset, H ):
    """
    Return the pols array for this step.
        - check is all the toltb has this pol.
    The order is:
    * local step value
    * global value
    * default = []
    """
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Pol" ] )
    # local val
    pols = parset.getStringVector( stepOptName, [] )
    if pols == []:
        pols = parset.getStringVector( "LoSoTo.Pols", [] )
        # global val or default

    return pols


def openSoltabs( H, ss_sts ):
    """
    Return a list of soltab objects checking the Soltab and SolType parameters
    """
    soltabs = []
    for ss_st in ss_sts:
        ss, st = ss_st.split('/')
        soltabs.append( H.getSoltab(ss, st) )

    return soltabs


def formatSelection(ant=[], pol=[], dir=[]):
    """
    return a string that can be used by the solFetcher as a selection criteria
    """
    s=''
    if ant != []:
        for a in ant:
            if a == ant[0]: s = s+'('
            else: s = s+' & '
            s = s+'(ant == \''+a+'\')'
            if a == ant[-1]: s = s+')'
    if pol != []:
        for p in pol:
            if p == pol[0]: s = s+'('
            else: s = s+' & '
            s = s+'(pol == \''+p+'\') '
            if p == pol[-1]: s = s+')'
    if dir != []:
        for d in dir:
            if d == dir[0]: s = s+'('
            else: s = s+' & '
            s = s+'(dir == \''+d+'\') '
            if d == dir[-1]: s = s+')'

    return s
