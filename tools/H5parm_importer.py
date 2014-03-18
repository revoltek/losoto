#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool is used to convert a parmdb into a H5parm format.
# It can be run on a globaldb created with parmdb_collector.py
# or on a single SB which contains the necessary
# sky/instrument/ANTENNA/FIELD tables.
# It handles Gain/DirectionalGain/RotationAngle/CommonRotationAngle/CommonScalarPhase solution types.

# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (fdg@hs.uni-hamburg.de)"

import sys, os, glob
import socket
import numpy as np
import progressbar
import logging
import pyrap.tables as pt
import lofar.parmdb
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm, solWriter


def parmdbToAxes(solEntry):
    """
    Extract the information written as a string in the parmdb format
    """
    pol = None; pol1 = None; pol2 = None;
    dir = None; ant = None; parm = None

    thisSolType = solEntry.split(':')[0]

    # For CommonRotationAngle assuming [CommonRotationAngle:ant]
    if thisSolType == 'CommonRotationAngle':
        thisSolType, ant = solEntry.split(':')
        dir = 'pointing'

    # For RotationAngle assuming [RotationAngle:ant:sou]
    elif thisSolType == 'RotationAngle':
        thisSolType, ant, dir = solEntry.split(':')

    # For TEC assuming [TEC:ant]
    elif thisSolType == 'TEC':
        thisSolType, ant = solEntry.split(':')
        dir = 'pointing'

    # For Clock assuming [Clock:ant]
    elif thisSolType == 'Clock':
        thisSolType, ant = solEntry.split(':')

    # For CommonScalarPhase assuming [CommonScalarPhase:ant]
    elif thisSolType == 'CommonScalarPhase':
        thisSolType, ant = solEntry.split(':')
        dir = 'pointing'

    # For ScalarPhase assuming [ScalarPhase:ant:sou]
    elif thisSolType == 'ScalarPhase':
        thisSolType, ant, dir = solEntry.split(':')

    # For Gain assuming [Gain:pol1:pol2:parm:ant]
    elif thisSolType == 'Gain':
        thisSolType, pol1, pol2, parm, ant = solEntry.split(':')
        dir = 'pointing'

    # For DirectionalGain assuming [DirecitonalGain:pol1:pol2:parm:ant:sou]
    elif thisSolType == 'DirectionalGain':
        thisSolType, pol1, pol2, parm, ant, dir = solEntry.split(':')

    else:
        logging.error('Unknown solution type "'+thisSolType+'". Ignored.')

    if pol1 != None and pol2 != None:
        if pol1 == '0' and pol2 == '0': pol = 'XX'
        if pol1 == '1' and pol2 == '0': pol = 'YX'
        if pol1 == '0' and pol2 == '1': pol = 'XY'
        if pol1 == '1' and pol2 == '1': pol = 'YY'

    return pol, dir, ant, parm


if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v] <H5parm> <globaldb/SBname> \n'\
                    +_author, version='%prog '+losoto._version.__version__)
    opt.add_option('-v', '--verbose', help='Go Vebose! (default=False)', action='store_true', default=False)
    opt.add_option('-s', '--solset', help='Solution-set name (default=sol###)', type='string', default=None)
    opt.add_option('-i', '--instrument', help='Name of the instrument table (default=instrument*)', type='string', default='instrument*')
    opt.add_option('-c', '--complevel', help='Compression level from 0 (no compression, fast) to 9 (max compression, slow) (default=5)', type='int', default='5')
    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()
    if options.verbose: losoto._logging.setLevel("debug")

    h5parmFile = args[0]
    logging.info("H5parm filename = "+h5parmFile)

    globaldbFile = args[1]
    if not os.path.exists(globaldbFile):
        logging.critical('Input globaldb/SB file not found.')
        sys.exit(1)
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

    # Make a list of all available instrument tables (only 1 for a standard MS)
    instrumentdbFiles = [ instrumentdbFile for instrumentdbFile in \
        glob.glob(os.path.join(globaldbFile,options.instrument)) \
        if os.path.isdir(instrumentdbFile) ]

    # open/create the h5parm file and the solution-set
    h5parm = h5parm(h5parmFile, readonly = False, complevel = complevel)

    solsetName = options.solset
    solset = h5parm.makeSolset(solsetName)

    # Create tables using the first instrumentdb
    # TODO: is there a better solution which check all the instrumentdbs?
    pdb = lofar.parmdb.parmdb(instrumentdbFiles[0])

    solTypes = list(set(x[0] for x in  (x.split(":") for x in pdb.getNames())))
    logging.info('Found solution types: '+', '.join(solTypes))

    # rewrite solTypes in order to put together
    # Gain <-> DirectionalGain
    # CommonRotationAngle <-> RotationAngle
    # CommonScalarPhase <-> ScalarPhase
    # it also separate Real/Imag/Ampl/Phase into different solTypes
    if "Gain" in solTypes:
        solTypes.remove('Gain')
        solTypes.append('*Gain:*:Real')
        solTypes.append('*Gain:*:Imag')
        solTypes.append('*Gain:*:Ampl')
        solTypes.append('*Gain:*:Phase')
    if "DirectionalGain" in solTypes:
        solTypes.remove('DirectionalGain')
        solTypes.append('*Gain:*:Real')
        solTypes.append('*Gain:*:Imag')
        solTypes.append('*Gain:*:Ampl')
        solTypes.append('*Gain:*:Phase')
    if "RotationAngle" in solTypes:
        solTypes.remove('RotationAngle')
        solTypes.append('*RotationAngle')
    if "CommonRotationAngle" in solTypes:
        solTypes.remove('CommonRotationAngle')
        solTypes.append('*RotationAngle')
    if "ScalarPhase" in solTypes:
        solTypes.remove('ScalarPhase')
        solTypes.append('*ScalarPhase')
    if "CommonScalarPhase" in solTypes:
        solTypes.remove('CommonScalarPhase')
        solTypes.append('*ScalarPhase')
    solTypes = list(set(solTypes))

    # every soltype creates a different solution-table
    for solType in solTypes:

        if len(pdb.getNames(solType+':*')) == 0: continue

        pols = set(); dirs = set(); ants = set();
        freqs = set(); times = set(); ptype = set()

        logging.info('Reading '+solType+'.')
        pbar = progressbar.ProgressBar(maxval=len(instrumentdbFiles)*len(pdb.getNames(solType+':*'))).start()
        ipbar = 0

        for instrumentdbFile in instrumentdbFiles:

            pdb = lofar.parmdb.parmdb(instrumentdbFile)

            # create the axes grid, necessary if not all entries have the same axes lenght
            data = pdb.getValuesGrid(solType+':*')
            for solEntry in data:

                pol, dir, ant, parm = parmdbToAxes(solEntry)

                if pol != None: pols |= set([pol])
                if dir != None: dirs |= set([dir])
                if ant != None: ants |= set([ant])
                freqs |= set(data[solEntry]['freqs'])
                times |= set(data[solEntry]['times'])
                pbar.update(ipbar)
                ipbar += 1

        pbar.finish()

        pols = np.sort(list(pols)); dirs = np.sort(list(dirs)); ants = np.sort(list(ants)); freqs = np.sort(list(freqs)); times = np.sort(list(times))
        shape = [i for i in (len(pols), len(dirs), len(ants), len(freqs), len(times)) if i != 0]
        vals = np.empty(shape)
        vals[:] = np.nan
        weights = np.zeros(shape)

        logging.info('Filling table.')
        pbar = progressbar.ProgressBar(maxval=len(instrumentdbFiles)*len(pdb.getNames(solType+':*'))).start()
        ipbar = 0

        for instrumentdbFile in instrumentdbFiles:

            pdb = lofar.parmdb.parmdb(instrumentdbFile)

            # fill the values
            data = pdb.getValuesGrid(solType+':*')
            if 'Real' in solType: dataIm = pdb.getValuesGrid(solType.replace('Real','Imag')+':*')
            if 'Imag' in solType: dataRe = pdb.getValuesGrid(solType.replace('Imag','Real')+':*')
            for solEntry in data:

                pol, dir, ant, parm = parmdbToAxes(solEntry)
                ptype |= set([solEntry.split(':')[0]]) # original parmdb solution type

                freq = data[solEntry]['freqs']
                time = data[solEntry]['times']

                val = data[solEntry]['values']

                # convert Real and Imag in Amp and Phase respectively
                if parm == 'Real':
                    solEntryIm = solEntry.replace('Real','Imag')
                    valI = dataIm[solEntryIm]['values']
                    val = np.sqrt((val**2)+(valI**2))
                if parm == 'Imag':
                    solEntryRe = solEntry.replace('Imag','Real')
                    valR = dataRe[solEntryRe]['values']
                    val = np.arctan2(val, valR)

                coords = []
                if pol != None:
                    polCoord = np.searchsorted(pols, pol)
                    coords.append(polCoord)
                if dir != None:
                    dirCoord = np.searchsorted(dirs, dir)
                    coords.append(dirCoord)
                if ant != None:
                    antCoord = np.searchsorted(ants, ant)
                    coords.append(antCoord)
                freqCoord = np.searchsorted(freqs, freq)
                timeCoord = np.searchsorted(times, time)
                vals[tuple(coords)][np.ix_(freqCoord,timeCoord)] = val.T
                weights[tuple(coords)][np.ix_(freqCoord,timeCoord)] = 1
                pbar.update(ipbar)
                ipbar += 1

        pbar.finish()
        if solType == '*RotationAngle':
            h5parm.makeSoltab(solset, 'rotation', axesNames=['dir','ant','freq','time'], \
                    axesVals=[dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == '*ScalarPhase':
            h5parm.makeSoltab(solset, 'scalarphase', axesNames=['dir','ant','freq','time'], \
                    axesVals=[dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == 'Clock':
            h5parm.makeSoltab(solset, 'clock', axesNames=['ant','freq','time'], \
                    axesVals=[ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == 'TEC':
            h5parm.makeSoltab(solset, 'tec', axesNames=['dir','ant','freq','time'], \
                    axesVals=[dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == '*Gain:*:Real' or solType == '*Gain:*:Ampl':
            h5parm.makeSoltab(solset, 'amplitude', axesNames=['pol','dir','ant','freq','time'], \
                    axesVals=[pols,dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))
        elif solType == '*Gain:*:Imag' or solType == '*Gain:*:Phase':
            h5parm.makeSoltab(solset, 'phase', axesNames=['pol','dir','ant','freq','time'], \
                    axesVals=[pols,dirs,ants,freqs,times], vals=vals, weights=weights, parmdbType=', '.join(list(ptype)))

    logging.info('Collecting information from the ANTENNA table.')
    antennaTable = pt.table(antennaFile, ack=False)
    antennaNames = antennaTable.getcol('NAME')
    antennaPositions = antennaTable.getcol('POSITION')
    antennaTable.close()
    antennaTable = solset._f_get_child('antenna')
    antennaTable.append(zip(*(antennaNames,antennaPositions)))

    logging.info('Collecting information from the FIELD table.')
    fieldTable = pt.table(fieldFile, ack=False)
    phaseDir = fieldTable.getcol('PHASE_DIR')
    pointing = phaseDir[0, 0, :]
    fieldTable.close()

    sourceTable = solset._f_get_child('source')
    # add the field centre, that is also the direction for Gain and CommonRotationAngle
    sourceTable.append([('pointing',pointing)])

    dirs = []
    for tab in solset._v_children:
        c = solset._f_getChild(tab)
        if c._v_name != 'antenna' and c._v_name != 'source':
            dirs.extend(list(set(c.dir)))
    # remove duplicates
    dirs = list(set(dirs))
    # remove any pointing (already in the table)
    if 'pointing' in dirs:
        dirs.remove('pointing')

    if dirs != []:
        logging.info('Collecting information from the sky table.')
        sourceFile = skydbFile + '/SOURCES'
        src_table = pt.table(sourceFile, ack=False)
        sub_tables = src_table.getsubtables()
        vals = []
        ra = dec = np.nan
        has_patches_subtable = False
        for sub_table in sub_tables:
            if 'PATCHES' in sub_table:
                has_patches_subtable = True
        if has_patches_subtable:
            # Read values from PATCHES subtable
            src_table.close()
            sourceFile = skydbFile + '/SOURCES/PATCHES'
            src_table = pt.table(sourceFile, ack=False)
            patch_names = src_table.getcol('PATCHNAME')
            patch_ras = src_table.getcol('RA')
            patch_decs = src_table.getcol('DEC')
            for source in dirs:
                try:
                    patch_indx = patch_names.index(source)
                    ra = patch_ras[patch_indx]
                    dec = patch_decs[patch_indx]
                except ValueError:
                    ra = np.nan
                    dec = np.nan
                    logging.error('Cannot find the source '+source+'. I leave NaNs.')
                vals.append([ra, dec])
            src_table.close()
        else:
            # Try to read default values from parmdb instead
            skydb = lofar.parmdb.parmdb(skydbFile)
            vals = []
            ra = dec = np.nan

            for source in dirs:
                try:
                    ra = skydb.getDefValues('Ra:' + source)['Ra:' + source][0][0]
                    dec = skydb.getDefValues('Dec:' + source)['Dec:' + source][0][0]
                except KeyError:
                    # Source not found in skymodel parmdb, try to find components
                    logging.warning('Cannot find the source '+source+'. Trying components.')
                    ra = np.array(skydb.getDefValues('Ra:*' + source + '*').values())
                    dec = np.array(skydb.getDefValues('Dec:*' + source + '*').values())
                    if len(ra) == 0 or len(dec) == 0:
                        ra = np.nan
                        dec = np.nan
                        logging.error('Cannot find the source '+source+'. I leave NaNs.')
                    else:
                        ra = ra.mean()
                        dec = dec.mean()
                        logging.info('Found average direction for '+source+' at ra:'+str(ra)+' - dec:'+str(dec))
                vals.append([ra, dec])
        sourceTable.append(zip(*(dirs,vals)))

    logging.info("Total file size: "+str(int(h5parm.H.get_filesize()/1024./1024.))+" M.")

    # Add CREATE entry to history and print summary of tables if verbose
    soltabs = h5parm.getSoltabs(solset=solset)
    for st in soltabs:
        sw = solWriter(soltabs[st])
        sw.addHistory('CREATE (by H5parm_importer.py from %s:%s/%s)' % (socket.gethostname(), os.path.abspath(''), globaldbFile))
    if options.verbose:
        logging.info(str(h5parm))

    del h5parm
    logging.info('Done.')
