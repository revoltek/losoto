#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool is used to convert a parmdb into a H5parm format.
# It can be run on a globaldb created with parmdb_collector.py
# or on a single SB which contains the necessary
# sky/instrument/ANTENNA/FIELD tables.
# It handles Gain/DirectionalGain/RotationAngle/CommonRotationAngle solution types.

# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)"

import sys, os, glob
import numpy as np
import progressbar
import logging
import lofar.parmdb
import pyrap.tables as pt
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm

# Options
import optparse
opt = optparse.OptionParser(usage='%prog [-v] [-p H5parm] [-g globaldb/SBname] \n'\
                +_author, version='%prog '+losoto._version.__version__)
opt.add_option('-v', '--verbose', help='Go VeRbOsE! (default=False)', action='store_true', default=False)
opt.add_option('-p', '--h5parm', help='H5parm output file (default=global.h5)', type='string', default='global.h5')
opt.add_option('-g', '--globaldb', help='Globaldb/MS name (default=globaldb)', type='string', default='globaldb')
opt.add_option('-s', '--solset', help='Solution-set name (default=sol###)', type='string', default='')
opt.add_option('-c', '--complevel', help='Compression level from 0 (no compression, fast) to 9 (max compression, slow) (default=9)', type='int', default='9')
(options, args) = opt.parse_args()

if options.verbose: losoto._logging.setVerbose()

h5parmFile = options.h5parm
logging.info("H5parm filename = "+h5parmFile)
globaldbFile = options.globaldb
logging.info("globaldb filename = "+globaldbFile)
complevel = options.complevel

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

# open/create the h5parm file and the solution-set
h5parm = h5parm(h5parmFile, readonly = False, complevel = complevel)

solsetName = options.solset
solset = h5parm.makeSolset(solsetName)
logging.info("Solution-set name = "+solset._v_name)

# Create tables using the first instrumentdb
# TODO: is there a better solution which check all the instrumentdbs?
pdb = lofar.parmdb.parmdb(instrumentdbFiles[0])

solTypes = set(x[0] for x in  (x.split(":") for x in pdb.getNames()))
logging.info('Found solution types: '+', '.join(solTypes))

# Fill the rotation table
if 'RotationAngle' in solTypes or 'CommonRotationAngle' in solTypes:
    soltabRot = h5parm.makeSoltab(solset, 'rotation', \
            descriptor=np.dtype([('time', np.float64),('freq',np.float64),('ant', np.str_, 16),('dir', np.str_, 16),('flag', np.bool),('val', np.float64)]))
    
    logging.info('Filling table...')
    pbar = progressbar.ProgressBar(maxval=len(instrumentdbFiles)*len(pdb.getNames('*RotationAngle:*'))).start()
    ipbar = 0
    for instrumentdbFile in instrumentdbFiles:
        pdb = lofar.parmdb.parmdb(instrumentdbFile)
        for solEntry in pdb.getNames('*RotationAngle:*'):
            solType = solEntry.split(':')[0]
    
            # For CommonRotationAngle assuming [CommonRotationAngle:ant]
            if solType == 'CommonRotationAngle':
                solType, ant = solEntry.split(':')
                val = pdb.getValuesGrid(solEntry)[solEntry]['values']
                parm = 'rotation'
                dir = 'pointing'

            # For RotationAngle assuming [RotationAngle:ant:sou]
            elif solType == 'RotationAngle':
                solType, ant, dir = solEntry.split(':')
                val = pdb.getValuesGrid(solEntry)[solEntry]['values']
                parm = 'rotation'

            else:
                logging.error('Unknown solution type "'+solType+'". Ignored.')

            # cycle on time and freq and add the value to the h5parm table
            times = pdb.getValuesGrid(solEntry)[solEntry]['times']
            freqs = pdb.getValuesGrid(solEntry)[solEntry]['freqs']
            pbar.update(ipbar)
            ipbar += 1
                
            # to speed up iterate only on freq
            for idxFreq, freq in enumerate(freqs):
                rows = zip(*(times, [freq]*len(times), [ant]*len(times), [dir]*len(times), [False]*len(times), val[:,idxFreq]))
                h5parm.addRow(soltabRot, rows)
    pbar.finish()

    # Index columns
    logging.info('Indexing columns...')
    for c in ['ant','freq','dir','time']:
        col = soltabRot.colinstances[c]
        col.create_index()

# fil amplitude and phase tables
if 'Gain' in solTypes or 'DirectionalGain' in solTypes:

    solParms = set(x[3] for x in  (x.split(":") for x in pdb.getNames('*Gain:*')))
    logging.info('Found Gain solution parameters: '+', '.join(solParms))

    if 'Ampl' in solParms or 'Imag' in solParms or 'Real' in solParms :
        soltabAmp = h5parm.makeSoltab(solset, 'amplitude', \
                descriptor=np.dtype([('time', np.float64),('freq',np.float64),('ant', np.str_, 16),('dir', np.str_, 16),('pol', np.str_, 2),('flag', np.bool),('val', np.float64)]))
    
    if 'Phase' in solParms or 'Imag' in solParms or 'Real' in solParms :
        soltabPhase = h5parm.makeSoltab(solset, 'phase', \
                descriptor=np.dtype([('time', np.float64),('freq',np.float64),('ant', np.str_, 16),('dir', np.str_, 16),('pol', np.str_, 2),('flag', np.bool),('val', np.float64)]))

    logging.info('Filling tables...')
    pbar = progressbar.ProgressBar(maxval=len(instrumentdbFiles)*len(pdb.getNames('*Gain:*'))).start()
    ipbar = 0
    for instrumentdbFile in instrumentdbFiles:
        pdb = lofar.parmdb.parmdb(instrumentdbFile)
        for solEntry in pdb.getNames('*Gain:*'):
            solType = solEntry.split(':')[0]
            # For Gain assuming [Gain:pol1:pol2:parm:ant]
            if solType == 'Gain':
                solType, pol1, pol2, parm, ant = solEntry.split(':')
                dir = 'pointing'

            # For DirectionalGain assuming [DirecitonalGain:pol1:pol2:parm:ant:sou]
            elif solType == 'DirectionalGain':
                solType, pol1, pol2, parm, ant, dir = solEntry.split(':')

            else:
                logging.error('Unknown solution type "'+solType+'". Ignored.')

            val = pdb.getValuesGrid(solEntry)[solEntry]['values']
            if pol1 == '0' and pol2 == '0': pol = 'XX'
            if pol1 == '1' and pol2 == '0': pol = 'YX'
            if pol1 == '0' and pol2 == '1': pol = 'XY'
            if pol1 == '1' and pol2 == '1': pol = 'YY'
            if parm == 'Ampl': parm = 'amplitude'
            if parm == 'Phase': parm = 'phase'
            # ugly workaround, when find a Real, retreive the amp and set it
            # when find an Imag, retreive the phase and set it
            if parm == 'Real':
                solEntry = solEntry.replace('Real','Imag')
                valI = pdb.getValuesGrid(solEntry)[solEntry]['values']
                val = np.sqrt((val**2)+(valI**2))
                parm = 'amplitude'
            if parm == 'Imag':
                solEntry = solEntry.replace('Imag','Real')
                valR = pdb.getValuesGrid(solEntry)[solEntry]['values']
                val = np.arctan2(val, valR)
                parm = 'phase'
        
            # cycle on time and freq and add the value to the h5parm table
            times = pdb.getValuesGrid(solEntry)[solEntry]['times']
            freqs = pdb.getValuesGrid(solEntry)[solEntry]['freqs']
            pbar.update(ipbar)
            ipbar += 1
            for idxFreq, freq in enumerate(freqs):
                rows = zip(*(times, [freq]*len(times), [ant]*len(times), [dir]*len(times), [pol]*len(times), [False]*len(times), val[:,idxFreq]))
                if parm == 'amplitude':
                    h5parm.addRow(soltabAmp, rows)
                else:
                    h5parm.addRow(soltabPhase, rows)

    pbar.finish()

    # Index columns
    logging.info('Indexing columns...')
    for c in ['ant','freq','pol','dir','time']:
        col = soltabAmp.colinstances[c]
        col.create_index()
        col = soltabPhase.colinstances[c]
        col.create_index()


logging.info('Collecting informations from the ANTENNA table.')
antennaTable = pt.table(antennaFile)
antennaNames = antennaTable.getcol('NAME')
antennaPositions = antennaTable.getcol('POSITION')
antennaTable.close()
descriptor = np.dtype([('name', np.str_, 16),('position', np.float32, 3)])
antennaTable = h5parm.H.createTable(solset, 'antenna', descriptor, 'Antenna names and positions')
antennaTable.append(zip(*(antennaNames,antennaPositions)))

logging.info('Collecting informations from the FIELD table.')
fieldTable = pt.table(fieldFile)
phaseDir = fieldTable.getcol('PHASE_DIR')
pointing = phaseDir[0, 0, :]
fieldTable.close()

descriptor = np.dtype([('name', np.str_, 16),('dir', np.float32, 2)])
sourceTable = h5parm.H.createTable(solset, 'source', descriptor, 'Source names and directions')
# add the field centre, that is also the direction for Gain and CommonRotationAngle
sourceTable.append([('pointing',pointing)])

dirs = []
for tab in solset._v_leaves:
    c = solset._f_getChild(tab)
    if c._v_name != 'antenna' and c._v_name != 'source':
        dirs.extend(c.col('dir'))
dirs = set(dirs)
dirs.remove('pointing')

if dirs != []:
    logging.info('Collecting informations from the sky table.')
    skydb = lofar.parmdb.parmdb(skydbFile)
    vals = []
    ra = dec = np.nan
    for source in set(dirs):
        try:
            ra = skydb.getDefValues('Ra:' + source)['Ra:' + source][0][0]
            dec = skydb.getDefValues('Dec:' + source)['Dec:' + source][0][0]
        except KeyError:
            # Source not found in skymodel parmdb, try to find components
            logging.error('Cannot find the source '+source+'. Trying components.')
            ra = np.array(skydb.getDefValues('Ra:*' + source + '*').values()).mean()
            dec = np.array(skydb.getDefValues('Dec:*' + source + '*').values()).mean()
            if ra == np.nan or dec == np.nan:
                logging.error('Cannot find the source '+source+'. I leave NaNs.')
            else:
                logging.info('Found average direction for '+source+' at ra:'+str(ra)+' - dec:'+str(dec))
        vals.append([ra, dec])
    sourceTable.append(zip(*(dirs,vals)))

logging.info("Total file size: "+str(h5parm.H.get_filesize()/1024./1024.)+" M")
del h5parm
logging.info('Done.')
