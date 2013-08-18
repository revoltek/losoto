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
import _version
import _logging
from h5parm import h5parm
import lofar.parameterset

# Options
import optparse
opt = optparse.OptionParser(usage='%prog [-v] h5parm parset [default: losoto.parset] \n'
        +_author, version='%prog '+_version.__version__)
opt.add_option('-v', '--verbose', help='Go VeRbOsE! (default=False)', action='store_true', default=False)
(options, args) = opt.parse_args()

if options.verbose: _logging.setVerbose()

# Check options
try: h5parmFile = args[0]
except:
    logging.critical('Missing H5parm file.')
    sys.exit(1)
try: parsetFile = args[1]
except: parsetFile = 'losoto.parset'

if not os.path.isfile(h5parmFile):
    logging.critical("Missing h5parm file.")
    sys.exit(1)
if not os.path.isfile(parsetFile):
    logging.critical("Missing parset file, I don't know what to do :'(")
    sys.exit(1)

# Open the H5parm
H = h5parm(h5parmFile, readonly=False)

# from ~vdtol/Expion-2011-05-03/src
parset = lofar.parameterset.parameterset( parsetFile )
steps = parset.getStringVector( "LoSoTo.Steps", [] )

# Possible operations, linked to relative function
import operations
operations = { "RESET": operations.reset,
               "EXAMPLE": operations.example,
#               "CLOCKTEC": operations.clocktec,
#               "FLAG": operations.flag,
               "SMOOTH": operations.smooth
#               "INTERP": operations.interp,
#               "PLOT": operations.plot,
#               "APPLY": operations.apply
}

for step in steps:
   operation = parset.getString( '.'.join( [ "LoSoTo.Steps", step, "Operation" ] ) )
   logging.info("--> Starting \'" + step + "\' step (operation: " + operation \
           + ").")
   returncode = operations[ operation ].run( step, parset, H )
   if returncode != 0:
      logging.error("Step \'" + step + "\' incomplete. Try to continue anyway.")
   else:
      logging.info("Step \'" + step + "\' completed successfully.")

del H
logging.info("Done.")
