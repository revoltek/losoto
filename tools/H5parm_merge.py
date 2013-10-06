#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool is used to merge multiple H5parm.

# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)"

import sys, os, glob
import numpy as np
import logging
import lofar.parmdb
import pyrap.tables as pt
import losoto._version
import losoto._logging
import losoto.h5parm

# Options
import optparse
opt = optparse.OptionParser(usage='%prog [-v] [-f H5parm:solset] [-t H5parm:solset] \n'\
                        +_author, version='%prog '+losoto._version.__version__)
opt.add_option('-v', '--verbose', help='Go VeRbOsE! (default=False)', action='store_true', default=False)
opt.add_option('-o', '--orig', help='H5parm origin file (filename:solset)', type='string', default=None)
opt.add_option('-d', '--dest', help='H5parm destination file (filename:solset)', type='string', default=None)
(options, args) = opt.parse_args()

if options.verbose: losoto._logging.setVerbose("info")

h5parmFrom = options.orig
h5parmTo = options.dest
if h5parmFrom == None or h5parmTo == None:
    logging.critical("Missing H5parm file.")
    sys.exit(1)

logging.info("H5parm origin = "+h5parmFrom)
logging.info("H5parm destination = "+h5parmTo)

# scompose input values
h5parmFromFile, solsetFrom = h5parmFrom.split(':')
h5parmToFile, solsetTo = h5parmTo.split(':')

# retrieve table
hf = losoto.h5parm.h5parm(h5parmFromFile)
ssF = hf.getSolset(solset=solsetFrom)

# write table
ht = losoto.h5parm.h5parm(h5parmToFile, readonly=False)
# check if the solset exists
if solsetTo in ht.getSolsets():
    logging.critical('Destination solset already exists, quitting.')
    sys.exit(1)
ssT = ht.makeSolset(solsetName = solsetTo)
# write the soltabs
ssF._f_copy_children(ssT)

del hf
del ht

logging.info("Done.")
