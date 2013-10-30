#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool is used to convert an H5parm file to parmdb format by writing to
# existing parmdb instrument table(s).
#
# It handles Gain/DirectionalGain/RotationAngle/CommonRotationAngle/CommonScalarPhase solution types.
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de), David Rafferty (drafferty@hs.uni-hamurg.de)"

import sys, os, glob, re
import numpy as np
import shutil
import progressbar
import logging
import lofar.parmdb
import pyrap.tables as pt
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm, solWriter, solFetcher


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

    # For CommonScalarPhase assuming [CommonScalarPhase:ant]
    elif thisSolType == 'CommonScalarPhase':
        thisSolType, ant = solEntry.split(':')

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

    if pol != None:
        pol = re.escape(pol)
    if dir != None:
        dir = re.escape(dir)
    if ant != None:
        ant = re.escape(ant)
    return pol, dir, ant, parm


def getSoltabFromSolType(solType, solTabs, parm='ampl'):
    """Return a list of solution tables that corresponds to the type given

    solType - string defining solution type. E.g., "DirectionalGain"
    solTabs - solution tables returned by h5parm.getSoltabs()
    parm - parm to return if solType is "Gain": 'ampl'/'real' or
           'phase'/'imag'

    The soltab parmdb_type attribute is used as an additional filter to
    to distinguish multiple possible matches. If it is not available, all
    matches are returned.
    """
    solTabList = []
    if parm != None:
        parm = parm.lower()

    for name, st in solTabs.iteritems():
        # Handle gain table separately, as we need to distinguish ampl and phase
        if solType == 'DirectionalGain' or solType == 'Gain':
            if (parm == 'ampl' or parm == 'real') and st._v_title == 'amplitude':
                if hasattr(st._v_attrs, 'parmdb_type'):
                    if solType in st._v_attrs['parmdb_type'].split(', '):
                        solTabList.append(st)
                else:
                    solTabList.append(st)
            elif (parm == 'phase' or parm == 'imag') and st._v_title == 'phase':
                if hasattr(st._v_attrs, 'parmdb_type'):
                    if solType in st._v_attrs['parmdb_type'].split(', '):
                        solTabList.append(st)
                else:
                    solTabList.append(st)
        else:
            if hasattr(st._v_attrs, 'parmdb_type'):
                if solType in st._v_attrs['parmdb_type'].split(', '):
                    solTabList.append(st)
            else:
                if (solType == 'RotationAngle' or solType == 'CommonRotationAngle') and st._v_title == 'rotation':
                    solTabList.append(st)
                elif solType == 'CommonScalarPhase' and st._v_title == 'scalarphase':
                    solTabList.append(st)

    if len(solTabList) == 0:
        return None
    else:
        return solTabList


if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog <H5parm filename> <input globaldb/SB filename>\n'+
        _author, version='%prog '+losoto._version.__version__)
    opt.add_option('-v', '--verbose', help='Go VeRbOsE!',
        action='store_true', default=False)
    opt.add_option('-s', '--solset', help='Name of solution set to export '
        '(default=sol000)', type='string', default='sol000')
    opt.add_option('-o', '--outfile', help='Filename of global/SB to export parmdb to '
        '(default=input globaldb/SB filename)', type='string', default=None)
    opt.add_option('-r', '--root', help='Root string to prepend to input parmdb instrument directories '
        '(default=solution-set name)', type='string', default=None)
    opt.add_option('-c', '--clobber', help='Clobber exising files '
        '(default=False)', action='store_true', default=False)
    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()
    if options.verbose: losoto._logging.setLevel("debug")

    # Check input H5parm file
    h5parmFile = args[0]
    if not os.path.exists(h5parmFile):
        logging.critical('Input H5parm file not found.')
        sys.exit(1)
    logging.info("Input H5parm filename = "+h5parmFile)

    # Open the h5parm file and get solution set names
    h5parm_in = h5parm(h5parmFile, readonly = True)
    solsetNames = h5parm_in.getSolsets()

    # Check input parmdb file
    globaldbFile = args[1]
    if not os.path.exists(globaldbFile):
        logging.critical('Input globaldb/SB file not found.')
        sys.exit(1)
    logging.info("Input globaldb/SB filename = "+globaldbFile)

    # Check input solution set name
    solsetName = options.solset
    if solsetName not in solsetNames:
        logging.critical('The solution set "'+solsetName+'" was not found in input H5parm file.')
        sys.exit(1)
    logging.info("Solution set name = "+solsetName)
    solset = h5parm_in.getSolset(solsetName)

    # Make output parmdb directory if needed
    out_globaldbFile = options.outfile
    if out_globaldbFile is None:
        out_globaldbFile = globaldbFile
    if not os.path.exists(out_globaldbFile):
        os.mkdir(out_globaldbFile)
    logging.info("Output globaldb/SB filename = "+out_globaldbFile)

    # Check output root
    outroot = options.root
    if outroot is None:
        outroot = solsetName

    # Make a list of all available instrument tables (only 1 for a standard MS)
    instrumentdbFiles = [ instrumentdbFile for instrumentdbFile in \
        glob.glob(os.path.join(globaldbFile,'instrument*')) \
        if os.path.isdir(instrumentdbFile) ]
    if len(instrumentdbFiles) == 0:
        logging.critical('No parmdb table(s) found in input globaldb/SB file.')
        sys.exit(1)
    instrumentdbFiles.sort()

    # Find solution table types using the first instrumentdb
    # TODO: is there a better solution which check all the instrumentdbs?
    pdb = lofar.parmdb.parmdb(instrumentdbFiles[0])
    solTypes = list(set(x[0] for x in  (x.split(":") for x in pdb.getNames())))
    logging.info('Found solution types: '+', '.join(solTypes))
    solTypes = list(set(solTypes))
    solTabs = h5parm_in.getSoltabs(solset)

    # For each solType, select appropriate solutions and construct
    # the dictionary to pass to pdb.addValues()
    len_sol = {}
    for solType in solTypes:
        len_sol[solType] = len(pdb.getNames(solType+':*'))

    for instrumentdbFile in instrumentdbFiles:
        out_instrumentdbFile = out_globaldbFile + '/' + outroot + '_' + instrumentdbFile.split('/')[-1]
        logging.info('Filling '+out_instrumentdbFile+':')

        # Remove existing instrumentdb (if clobber) and create new one
        if os.path.exists(out_instrumentdbFile):
            if options.clobber:
                shutil.rmtree(out_instrumentdbFile)
            else:
                logging.critical('Output instrumentdb file exists and '
                    'clobber = False.')
                sys.exit(1)
        pdb_out = lofar.parmdb.parmdb(out_instrumentdbFile+'/', create=True)
        pbar = progressbar.ProgressBar(maxval=sum(len_sol.values())).start()
        ipbar = 0

        pdb_in = lofar.parmdb.parmdb(instrumentdbFile)
        for solType in solTypes:
            if len_sol[solType] == 0: continue

            solEntries = pdb_in.getNames(solType+':*')
            data = pdb_in.getValuesGrid(solType+':*')
            data_out = data.copy()
            for solEntry in solEntries:

                pol, dir, ant, parm = parmdbToAxes(solEntry)
                solTabList = getSoltabFromSolType(solType, solTabs, parm=parm)
                if solTabList is None:
                    logging.critical('Mismatch between parmdb table and H5parm '
                        'solution table: No solution tables found that match '
                        'the parmdb entry "'+solType+'"')
                    sys.exit(1)
                if len(solTabList) > 1:
                    logging.warning('More than one solution table found in H5parm '
                        'matching parmdb entry "'+solType+'". Taking the first match.')
                solTab = solTabList[0]
                sf = solFetcher(solTab)

                if pol == None and dir == None:
                    sf.setSelection(ant=ant)
                elif pol == None and dir != None:
                    sf.setSelection(ant=ant, dir=dir)
                elif pol != None and dir == None:
                    sf.setSelection(ant=ant, pol=pol)
                else:
                    sf.setSelection(ant=ant, pol=pol, dir=dir)

                # If needed, convert Amp and Phase to Real and Imag respectively
                if parm == 'Real':
                    SolTabList = getSoltabFromSolType(solType, solTabs, parm='phase')
                    soltab_phase = SolTabList[0]
                    sf_phase = solFetcher(soltab_phase, ant=ant, pol=pol, dir=dir)
                    val_amp = sf.getValues()[0]
                    val_phase = sf_phase.getValues()[0]
                    val = val_amp * np.cos(val_phase)
                elif parm == 'Imag':
                    SolTabList = getSoltabFromSolType(solType, solTabs, parm='ampl')
                    soltab_amp = SolTabList[0]
                    sf_amp = solFetcher(soltab_amp, ant=ant, pol=pol, dir=dir)
                    val_phase = sf.getValues()[0]
                    val_amp = sf_amp.getValues()[0]
                    val = val_amp * np.sin(val_phase)
                else:
                    val = sf.getValues()[0]

                # Match the frequency or frequencies of instrumentdb under
                # consideration
                sffreqs = sf.freq
                freqs = data[solEntry]['freqs']
                freq_list = [freq for freq in freqs if freq in sffreqs]
                shape = data_out[solEntry]['values'].shape
                if len(freq_list) == 0:
                    for i in range(len(val.shape)):
                        freq_ind_list.append(slice(None))
                    freq_ind = tuple(freq_ind_list)
                else:
                    freqAxisIdx = sf.getAxesNames().index('freq')
                    freq_ind = []
                    for i in range(len(val.shape)):
                        freq_ind.append(slice(None))
                    freq_ind[freqAxisIdx] = np.where(sffreqs == freq_list)
                    freq_ind = tuple(freq_ind)

                try:
                    data_out[solEntry]['values'] = val[freq_ind].reshape(shape)
                except ValueError, err:
                    logging.critical('Mismatch between parmdb table and H5parm '
                    'solution table: Differing number of frequencies and/or times')
                    sys.exit(1)
                pbar.update(ipbar)
                ipbar += 1
            pdb_out.addValues(data_out)

        pbar.finish()

    logging.info('Done.')






