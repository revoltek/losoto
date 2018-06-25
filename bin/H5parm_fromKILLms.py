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

    inputFile = args[1:]
    logging.info("KILLMS filenames = "+str(inputFile))
    h5parmFile = args[0]
    logging.info("H5parm filename = "+h5parmFile)
    
    # Common options
    complevel = options.complevel
    solsetName = options.solset

    for FileName in inputFile:

        SolsDico = np.load(FileName)
        Sols = SolsDico["Sols"]
        Sols = Sols.view(np.recarray)
        nt,nch,na,nd,_,_ = Sols.G.shape

        # build direction subtable
        ClusterCat = SolsDico["ClusterCat"]
        dirs = []
        for i, c in enumerate(ClusterCat):
            print type(c)
            dirs.append['Dir%02i' % i, c[1], c[2]]

        # build antenna subtable
        StationNames = SolsDico["StationNames"]


    print SolsDico.keys()
    print SolsDico['FreqDomains']
    print nt,nch,na,nd
    print tpye(Sols)
    print dirs
