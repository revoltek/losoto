#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This tool is used to merge 2 H5parms.
# Each h5parm solset will be copied into the last one.

_author = "Francesco de Gasperin (astro@voo.it)"

import sys, os, glob
import numpy as np
import pyrap.tables as pt
from losoto import _version
from losoto import _logging
import losoto.h5parm

if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v] <H5parm:solset> <H5parm:solset> \n'\
                            +_author, version='%prog '+_version.__version__)
    opt.add_option('-V', '--verbose', help='Go VERBOSE! (default=False)', action='store_true', default=False)
    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()

    # log
    logger = _logging.Logger('info')
    logging = _logging.logger
    if options.verbose: logger.set_level("debug")

    h5parmFrom = args[0]
    h5parmTo = args[1]

    logging.info("H5parm origin = "+h5parmFrom)
    logging.info("H5parm destination = "+h5parmTo)

    # scompose input values
    h5parmFromFile, solsetFrom = h5parmFrom.split(':')
    h5parmToFile, solsetTo = h5parmTo.split(':')

    if not (os.path.isfile(h5parmFromFile) and os.path.isfile(h5parmToFile)):
        logging.critical("Missing H5parm file.")
        sys.exit(1)

    # retrieve table
    hf = losoto.h5parm.h5parm(h5parmFromFile)
    ssF = hf.getSolset(solset=solsetFrom)

    # write table
    ht = losoto.h5parm.h5parm(h5parmToFile, readonly=False)
    # check if the solset exists
    if solsetTo in ht.getSolsetNames():
        logging.critical('Destination solset already exists, quitting.')
        sys.exit(1)
    ssT = ht.makeSolset(solsetName = solsetTo, addTables=False)
    # write the soltabs
    ssF.obj._f_copy_children(ssT.obj, recursive=True)

    del hf

    # Add entry to history
    soltabs = ht.getSolset(solsetTo).getSoltabs()
    for st in soltabs:
        st.addHistory('Copied (from {0}:{1})'.format(h5parmFromFile, solsetFrom))
    del ht

    logging.info("Done.")
