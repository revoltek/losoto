#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Authors:
# Francesco de Gasperin
# David Raffery
# Cyril Tasse
# Reinout van Weeren
# Maaijke Mevius
# Bas van der Tol
_author = "Francesco de Gasperin (fdg@hs.uni-hamburg.de)"

import sys
import os
import time
import atexit
import tables
import logging
import _version
import _logging
from h5parm import h5parm
import lofar.parameterset

def my_close_open_files(verbose):
    open_files = tables.file._open_files
    are_open_files = len(open_files) > 0
    if verbose and are_open_files:
        print >> sys.stderr, "Closing remaining open files:",
    # Compatibility fix
    if tables.__version__>='3.1.0':
        for fileh in list(open_files.handlers):
            if verbose:
                print >> sys.stderr, "%s..." % (fileh.filename,),
            fileh.close()
    else:
        for fileh in open_files.keys():
            if verbose:
                print >> sys.stderr, "%s..." % (open_files[fileh].filename,),
            open_files[fileh].close()
        if verbose:
            print >> sys.stderr, "done",
    if verbose and are_open_files:
        print >> sys.stderr

if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v|-V] h5parm parset [default: losoto.parset] \n'
            +_author, version='%prog '+_version.__version__)
    opt.add_option('-q', help='Quiet', action='store_true', default=False)
    opt.add_option('-v', help='Verbose', action='store_true', default=False)
    opt.add_option('-f', '--filter', help='Filter to use with "-i" option to filter on solution set names (default=None)', type='string', default=None)
    opt.add_option('-i', help='List information about h5parm file (default=False). A filter on the solution set names can be specified with the "-f" option.', action='store_true', default=False)
    (options, args) = opt.parse_args()

    atexit.register(my_close_open_files, False) # Suppress info about closing open files at exit
    if options.q:
        _logging.setLevel('warning')
    if options.v:
        _logging.setLevel('debug')
        atexit.register(my_close_open_files, True) # Print info about closing open files at exit

    # Check options
    if len(args) not in [1, 2]:
        opt.print_help()
        sys.exit()

    try: h5parmFile = args[0]
    except:
        logging.critical('Missing H5parm file.')
        sys.exit(1)

    try: parsetFile = args[1]
    except: parsetFile = 'losoto.parset'

    if not os.path.isfile(h5parmFile):
        logging.critical("Missing h5parm file.")
        sys.exit(1)
    if not os.path.isfile(parsetFile) and not options.i:
        logging.critical("Missing parset file, I don't know what to do :'(")
        sys.exit(1)

    # Open the H5parm
    H = h5parm(h5parmFile, readonly=False)

    # List h5parm information if desired
    if options.i:
        print(H.printInfo(options.filter))
        sys.exit()

    # from ~vdtol/Expion-2011-05-03/src
    parset = lofar.parameterset.parameterset( parsetFile )
    steps = parset.getStringVector( "LoSoTo.Steps", [] )

    # Possible operations, linked to relative function
    import operations
    operations = { 
                   "ABS": operations.abs,
                   "CLIP": operations.clip,
                   "INTERP": operations.interp,
                   "NORM": operations.norm,
                   "PLOT": operations.plot,
                   "RESET": operations.reset,
                   "SMOOTH": operations.smooth,
                   "TECFIT": operations.tecfit,
                   "TECSCREEN": operations.tecscreen,
                   # example operation
                   "EXAMPLE": operations.example
                   # TBI operations
                   #"CLOCKTEC": operations.clocktec,
                   #"FLAG": operations.flag,
                   #"COPY": operations.copy,
    }

    for step in steps:
       operation = parset.getString( '.'.join( [ "LoSoTo.Steps", step, "Operation" ] ) )
       logging.info("--> Starting \'" + step + "\' step (operation: " + operation + ").")
       start = time.clock()
       returncode = operations[ operation ].run( step, parset, H )
       if returncode != 0:
          logging.error("Step \'" + step + "\' incomplete. Try to continue anyway.")
       else:
          logging.info("Step \'" + step + "\' completed successfully.")
       elapsed = (time.clock() - start)
       logging.debug("Time for this step: "+str(elapsed)+" s.")


    del H
    logging.info("Done.")

