#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool is used to convert a parmdb into a H5parm format.
# It can be run on a globaldb created with parmdb_collector.py
# or on a single SB which contains the necessary
# sky/instrument/ANTENNA/FIELD tables.

# Authors:
# Francesco de Gasperin
# Bas van der Tol
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)"

import sys, os, glob
import numpy as np
import logging
import lofar.parmdb
import pyrap.tables as pt
import losoto._version
import losoto._logging

# Options
import optparse
opt = optparse.OptionParser(usage='%prog [-v] [-h H5parm] [-g globaldb/SBname] \n'\
                +_author, version='%prog '+losoto._version.__version__)
opt.add_option('-v', '--verbose', help='Go VeRbOsE! (default=False)', action='store_true', default=False)
opt.add_option('-p', '--h5parm', help='H5parm output file (default=global.h5)', type='string', default='global.h5')
opt.add_option('-g', '--globaldb', help='Output globaldb name (default=globaldb)', type='string', default='globaldb')
(options, args) = opt.parse_args()

if options.verbose: losoto._logging.setVerbose()

h5parmFile = options.h5parm
logging.info("H5parm filename = "+h5parmFile)
globaldbFile = options.globaldb
logging.info("globaldb filename = "+globaldbFile)

# Check is all the necessary files are available
antennaFile = os.path.join(globaldbFile,'ANTENNA')
if not os.path.isdir(antennaFile):
    logging.critical('Missing ANTENNA table.')
    sys.exit(1)
fieldFile = os.path.join(globaldbFile,'FIELD')
if not os.path.isdir(fieldFile):
    logging.critical('Missing FIELD table.')
    sys.exit(1)
skydbFile = os.path.join(globaldbFile,'sky')
if not os.path.isdir(skydbFile):
    logging.critical('Missing skydb table.')
    sys.exit(1)

# Make a list of all available instrument tables
instrumentdbFiles = [ instrumentdbFile for instrumentdbFile in \
    glob.glob(os.path.join(globaldbFile,'instrument*')) \
    if os.path.isdir(instrumentdbFile) ]

# Collect all the parameters from all the instrument tables
for instrumentdbFile in instrumentdbFiles:
    db = lofar.parmdb.parmdb(instrumentdbFile)
    for solEntry in db.getNames():
        solType = solEntry.split(':')[0]
        if solType == 'Gain':
            solType,pol1,pol2,parm,ant = np.array([ i.split(':') for i in db.getNames() ]).transpose()

        elif solType == 'DirectionalGain':

        elif solType == 'CommonRotationAngle':

        elif solType == 'RotationAngle':

        else:
            logging.error('Unknown solution type '+solType+'. Ignored.')
        
    gtype = list(set(gtypes))
    corrlist = list(set(corrs))
    parmlist = list(set(parms))
    antlist = list(set(ants))

print gtype, corrlist, parmlist, antlist


logging.info('Collecting informations from the ANTENNA table.')
antennaTable = pt.table(antennaFile)
antennaNames = antennaTable.getcol('NAME')
antennaPositions = antennaTable.getcol('POSITION')
antenna_table.close()

logging.info('Collecting informations from the FIELD table.')
fieldTable = pt.table(fieldFile)
phaseDir = fieldTable.getcol('PHASE_DIR')
pointing = phaseDir[0, 0, :]
fieldTable.close()

logging.info('Collecting informations from the sky table.')
skydb = lofar.parmdb.parmdb(skydbFile)
if self.DirectionalGainEnable:
            if len(sources) == 0:
                sources = ['*']
            self.sources = get_source_list(instrumentdb, sources)
            self.source_positions = []
            for source in self.sources:
                try:
                    RA = skydb.getDefValues('Ra:' + source)['Ra:' + source][0][0]
                    dec = skydb.getDefValues('Dec:' + source)['Dec:' + source][0][0]
                except KeyError:

                 # Source not found in skymodel parmdb, try to find components

                    RA = numpy.array(skydb.getDefValues('Ra:' + source + '.*').values()).mean()
                    dec = numpy.array(skydb.getDefValues('Dec:' + source + '.*').values()).mean()
                self.source_positions.append([RA, dec])
            else:
                self.sources = ['Pointing']
                self.source_positions = [list(self.pointing)]
            self.N_sources = len(self.sources)



