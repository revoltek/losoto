#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This tool is used to convert a parmdb into a H5parm format.
It can be run on a globaldb created with parmdb_collector.py
or on a single SB which contains the necessary
sky/instrument/ANTENNA/FIELD tables.
It can also load a dictionary from a json or pickle file with the keys 
instrumentdb (array with file names), antenna, field and skydb containing the 
corresponding file names.
It handles Gain/DirectionalGain/RotationAngle/
           Clock/TEC/CommonRotationAngle/CommonScalarPhase/CommonScalarAmpitude solution types.
"""
# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (astro@voo.it)"

import sys, os, glob
import socket
from losoto import _version
from losoto import _logging
from losoto._importer import create_h5parm

    
if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v] <H5parm> <globaldb/SBname> \n'\
                    +_author, version='%prog '+_version.__version__)
    opt.add_option('-V', '--verbose', help='Go Vebose! (default=False)', action='store_true', default=False)
    opt.add_option('-s', '--solset', help='Solution-set name (default=sol###)', type='string', default=None)
    opt.add_option('-i', '--instrument', help='Name of the instrument table (default=instrument*)', type='string', default='instrument*')
    opt.add_option('-c', '--complevel', help='Compression level from 0 (no compression, fast) to 9 (max compression, slow) (default=5)', type='int', default='5')
    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()

    # log
    logger = _logging.Logger('info')
    logging = _logging.logger
    if options.verbose: logger.set_level("debug")

    h5parmFile = args[0]
    logging.info("H5parm filename = "+h5parmFile)
    
    # Common options
    complevel = options.complevel
    solsetName = options.solset

    input_file = args[1]
    
    if input_file.endswith(".json"):
        try:
            import json
            logging.info("Loading json file: {}".format(input_file))
            files = json.load(open(input_file,"r"))
            instrumentdbFiles = [str(f) for f in files["instrumentdb"]]
            antennaFile = str(files["antenna"])
            fieldFile = str(files["field"])
            skydbFile = str(files["skydb"])
            globaldbFile = None
        except:
            logging.critical('Loading failed')
    elif input_file.endswith(".pckl"):
        try:
            import pickle
            logging.info("Loading pickle file: {}".format(input_file))
            files = pickle.load(open(input_file,"r"))
            instrumentdbFiles = files["instrumentdb"]
            antennaFile = files["antenna"]
            fieldFile = files["field"]
            skydbFile = files["skydb"]
            globaldbFile = None
        except:
            logging.critical('Loading failed')
    else:
        globaldbFile = input_file
        if not os.path.exists(globaldbFile):
            logging.critical('Input globaldb/SB file not found.')
            sys.exit(1)
        logging.info("globaldb filename = "+globaldbFile)
        
        # Load the path of the files
        antennaFile = os.path.join(globaldbFile,'ANTENNA')
        fieldFile = os.path.join(globaldbFile,'FIELD')
        skydbFile = os.path.join(globaldbFile,'sky')

        # Make a list of all available instrument tables (only 1 for a standard MS)
        instrumentdbFiles = [ instrumentdbFile for instrumentdbFile in \
            glob.glob(os.path.join(globaldbFile,options.instrument)) \
            if os.path.isdir(instrumentdbFile) ]

    # Check antennaFile and fieldFile
    if not os.path.isdir(antennaFile):
        logging.critical('Missing ANTENNA table.')
        sys.exit(1)
    if not os.path.isdir(fieldFile):
        logging.critical('Missing FIELD table.')
        sys.exit(1)
    
    # Call the method that creates the h5parm file
    create_h5parm(instrumentdbFiles, antennaFile, fieldFile, skydbFile,
                  h5parmFile, complevel, solsetName, globaldbFile=globaldbFile,verbose=options.verbose)
