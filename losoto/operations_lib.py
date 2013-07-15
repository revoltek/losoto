#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Some utilities for operations.py

import sys

# TODO: add pattern matching
def getStations( step, parset, H ):
    """
    Return the station array for this step.
    The order is:
    * local step value
    * global value
    * default
    """
    allStations = H.stations.names
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Stations" ] )
    # local val
    stations = parset.getStringVector( stepOptName, [] )
    if stations == []:
        # global val or default
        stations = parset.getStringVector( "LoSoTo.Stations", allStations )

    # sanity check
    for station in stations:
        if station not in allStations:
            print "ERROR: cannot find station", station, ", ignoring."
            stations.remove(station)
    return stations


def getSolTypes( step, parset, H ):
    """
    Return the solTypes array for this step.
    The order is:
    * local step value
    * global value
    * default
    """
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "SolType" ] )
    # local val
    solTypes = parset.getStringVector( stepOptName, [] )
    if solTypes == []:
        # global value or default
        solTypes = parset.getStringVector( "LoSoTo.SolTypes", ['Amp'] )

    # sanity check
    for solType in solTypes:
        if solType not in H.root:
            print "ERROR: solution type", solType, " in not in the HDF5 file, ignoring"
            solTypes.remove(solType)

    return solTypes


def getPols( step, parset, H ):
    """
    Return the pols array for this step.
    The order is:
    * local step value
    * global value
    * default
    """
    allPols = H.polarizations
    stepOptName = '.'.join( [ "LoSoTo.Steps", step, "Pols" ] )
    # local val
    pols = parset.getStringVector( stepOptName, [] )
    if pols == []:
        pols = parset.getStringVector( "LoSoTo.Pols", ['XX','YY'] )
        # global val or default

    # sanity check
    for pol in pols:
        if pol not in allPols:
            print "ERROR: unknown solution typr", solType, ", ignoring."
            pols.remove(pol)

    return pols
