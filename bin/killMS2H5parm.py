#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This tool is used to import killms solutions into a H5parm format.
"""
# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (astro@voo.it)"

import sys, os, glob, pickle
import logging
import numpy as np
from losoto import _version
from losoto import _logging
from losoto import h5parm as h5parm_mod
    
if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v] <H5parm> <killmsfile> \n'\
                    +_author, version='%prog '+_version.__version__)
    opt.add_option('-V', '--verbose', help='Go Vebose! (default=False)', action='store_true', default=False)
    opt.add_option('-s', '--solset', help='Solution-set name (default=sol###)', type='string', default=None)
    opt.add_option('-c', '--complevel', help='Compression level from 0 (no compression, fast) to 9 (max compression, slow) (default=5)', type='int', default='5')
    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()
    if options.verbose: _logging.setLevel("debug")

    inputFile = args[1]
    logging.info("KILLMS filenames = "+str(inputFile))
    h5parmFile = args[0]
    logging.info("H5parm filename = "+h5parmFile)
    
    # Common options
    complevel = options.complevel
    solsetName = options.solset

    SolsDico = np.load(inputFile)
    Sols = SolsDico["SolsTEC"]
    Sols = Sols.view(np.recarray)

    # build direction subtable
    ClusterCat = SolsDico["ClusterCat"]
    dirCoords = []; dirNames = []
    for i, c in enumerate(ClusterCat):
        dirCoords.append([c[1], c[2]])
        dirNames.append('Dir%02i' % i)

    # build antenna subtable
    stationNames = SolsDico["StationNames"]
    antPos = []; antNames = []
    for i, a in enumerate(stationNames):
        antPos.append([a,0,0,0])
        antNames.append(a)

    print SolsDico.keys()
    print Sols.dtype.names
    pols = ['XX','XY','YX','YY']
    times = (Sols["t0"]+Sols["t1"])/2.
    freqs = (SolsDico['FreqDomains'][0]+SolsDico['FreqDomains'][1])/2.
    vals = Sols['G']
    weights = np.ones(shape=vals.shape)
    print Sols['G'].shape

    h5parm = h5parm_mod.h5parm(h5parmFile, readonly = False, complevel = complevel)
    solset = h5parm.makeSolset(solsetName)
    solset.makeSoltab('phase', axesNames=['pol','dir','ant','freq','time'], \
            axesVals=[pols,dirNames,antNames,freqs,times], vals=vals, weights=weights)

    # fill source table
    sourceTable = solset.obj._f_get_child('source')
    sourceTable.append(zip(*(dirNames,dirCoords)))

    # fill antenna table
    antennaTable = solset.obj._f_get_child('antenna')
    antennaTable.append(zip(*(antNames,antPos)))






