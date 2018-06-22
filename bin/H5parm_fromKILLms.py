#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This tool is used to import killms solutions into a H5parm format.
"""
# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (astro@voo.it)"

import sys, os, glob
import logging
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

    input_file = args[1]
    h5parm_file = args[0]
    logging.info("H5parm filename = "+h5parmFile)
    
    # Common options
    complevel = options.complevel
    solset_name = options.solset


    
