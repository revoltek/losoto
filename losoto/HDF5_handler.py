#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Retrieving ang writing back data in H5parm format

def getSols( H, solType, station=[], pol=[], source=[] ):
    """
    Return solution array once selected for solType, station, pol and source.
    solType must be specified.
    """

    if solType == 'amp':

    elif solType == 'phase':

    elif solType == 'clock':

    elif solType == 'TEC':

    else:
        print "ERROR: unknown solType", solType, "."
        sys.exit(1)

