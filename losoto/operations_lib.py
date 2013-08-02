#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Some utilities for operations.py

import sys
import logging

# TODO: add pattern matching
def getParAnts( step, parset, H ):
    """
    Return the Ant array for this step.
    The order is:
    * local step value
    * global value
    * default = all
    """
    allAnts = H.getAnt(getParSolsets(step, parset, H)).keys()
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Ant" ] )
    # local val
    Ants = parset.getStringVector( stepOptName, [] )
    if Ants == []:
        # global val or default
        Ants = parset.getStringVector( "LoSoTo.Ant", allAnts )

    # sanity check
    for Ant in Ants:
        if Ant not in allAnts:
            logging.error("Cannot find Ant", Ant, ", ignoring.")
            Ants.remove(Ant)
    return Ants


def getParSoltabs( step, parset, H ):
    """
    Return the solution-table array for this step.
    The order is:
    * local step value
    * global value
    * default = all
    """
    allSoltabs = H.getSoltabs(getParSolsets(step, parset, H)).keys()
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Soltab" ] )
    # local val
    soltabs = parset.getStringVector( stepOptName, [] )
    if soltabs == []:
        # global value or default
        soltabs = parset.getStringVector( "LoSoTo.Soltab", allSoltabs )

    # sanity check
    for soltab in soltabs:
        if soltab not in allSoltabs:
            logging.error("Solution-table", soltab, " in not in the HDF5 file, ignoring")
            soltabs.remove(soltab)

    return soltabs


def getParSolsets( step, parset, H ):
    """
    Return the solution-set array for this step.
    The order is:
    * local step value
    * global value
    * default = all
    """
    allSolsets = H.getSolsets().keys()
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Solset" ] )
    # local val
    solsets = parset.getStringVector( stepOptName, [] )
    if solsets == []:
        # global value or default
        solsets = parset.getStringVector( "LoSoTo.Solset", allSolsets )

    # sanity check
    for solset in solsets:
        if solset not in allSolsets:
            logging.error("Solution-set", solset, " in not in the HDF5 file, ignoring")
            solsets.remove(solset)

    return solsets


def getParDirs( step, parset, H ):
    """
    Return the directions array for this step.
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


def openSoltabs( H, solsets, soltabs ):
    """
    Return a list of soltab objects
    """
    allSoltabs = []
    for solset in solsets:
        for soltab in soltabs:
            allSoltabs.append( H.getSoltab(solset, soltab) )

    return allSoltabs


def formatSelection(ant=[], pol=[], dir=[]):
    """
    return a string that can be used by the solFetcher as a selection criteria
    """
    s=''
    if ant != []:
        for a in ant:
            if a == ant[0]: s = s+'('
            else: s = s+' & '
            s = s+'(ant == '+a+')'
            if a == ant[-1]: s = s+')'
    if pol != []:
        for p in pol:
            if p == pol[0]: s = s+'('
            else: s = s+' & '
            s = s+'(pol == '+p+') '
            if p == pol[-1]: s = s+')'
    if ant != []:
        for d in dir:
            if d == dir[0]: s = s+'('
            else: s = s+' & '
            s = s+'(dir == '+d+') '
            if d == dir[-1]: s = s+')'

    return s












