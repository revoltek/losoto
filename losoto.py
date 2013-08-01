#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Authors:
# Francesco de Gasperin
# Cyril Tasse
# Reinout van Weeren
# Maaijke Mevius
# Bas van der Tol
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)"

import sys, os
import logging
import losoto.operations
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm

# Options
import optparse
opt = optparse.OptionParser(usage='%prog [-v] h5parm parset [default: losoto.parset] \n'
        +_author, version='%prog '+losoto._version.__version__)
opt.add_option('-v', '--verbose', help='Go VeRbOsE! (default=False)', action='store_true', default=False)
(options, args) = opt.parse_args()

if options.verbose: losoto._logging.setVerbose()

# Check options
try: parsetFile = args[0]
except: parsetFile = 'losoto.parset'
try: h5parmFile = args[1]
except:
    logging.critical('Missing H5parm file.')
    sys.exit(1)

if not os.path.isfile(parsetFile):
    logging.critical("Missing parset file, I don't know what to do :'(")
    sys.exit(1)

# Open the H5parm
H = h5parm(h5parmFile)

# from ~vdtol/Expion-2011-05-03/src
parset = lofar.parameterset.parameterset( parsetFile )
steps = parset.getStringVector( "LoSoTo.Steps", [] )

# Possible operations, linked to relative function
operations = { "RESET": losoto.operations.reset,
               "CLOCKTEC": losoto.operations.clocktec,
               "FLAG": losoto.operations.flag,
               "SMOOTH": losoto.operations.smooth,
               "INTERP": losoto.operations.interp,
               "PLOT": losoto.operations.plot,
               "APPLY": losoto.operations.apply }

for step in steps:
   operation = parset.getString( '.'.join( [ "LoSoTo.Steps", step, "Operation" ] ) )
   logging.info("--> Starting \'" + step + "\' step (operation: " + operation \
           + ") on solution-set "+solset+".")
   solset = operations[ operation ] ( step, parset, H, solset )

logging.info("Done.")
