#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This tool is used to convert an H5parm file to parmdb format by writing to
# existing parmdb instrument table(s).
#
# It handles Gain/DirectionalGain/RotationAngle/CommonRotationAngle/CommonScalarPhase solution types.
_author = "Francesco de Gasperin (astro@voo.it), David Rafferty (drafferty@hs.uni-hamburg.de)"

import sys, os, glob, re, time
import numpy as np
import shutil
import pyrap.tables as pt
import lofar.parmdb
from losoto import _version
from losoto import _logging
from losoto.h5parm import h5parm
try:
    import progressbar
except ImportError:
    import losoto.progressbar as progressbar


def parmdbToAxes(solEntry):
    """
    Extract the information written as a string in the parmdb format
    """
    pol = None
    pol1 = None
    pol2 = None
    dir = None
    ant = None
    parm = None

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
        dir = 'pointing'

    # For CommonScalarPhase assuming [CommonScalarPhase:ant]
    elif thisSolType == 'CommonScalarAmplitude':
        thisSolType, ant = solEntry.split(':')
        dir = 'pointing'

    # For ScalarPhase assuming [ScalarPhase:ant:sou]
    elif thisSolType == 'ScalarPhase':
        thisSolType, ant, dir = solEntry.split(':')

    # For ScalarPhase assuming [ScalarPhase:ant:sou]
    elif thisSolType == 'ScalarAmplitude':
        thisSolType, ant, dir = solEntry.split(':')

    # For TEC assuming [TEC:ant or TEC:pol:ant]
    elif thisSolType == 'TEC':
        try:
            thisSolType, ant = solEntry.split(':')
        except:
            thisSolType, pol, ant = solEntry.split(':')
            pol1 = pol
            pol2 = pol
        dir = 'pointing'

    # For Clock assuming [Clock:ant or Clock:pol:ant]
    elif thisSolType == 'Clock':
        try:
            thisSolType, ant = solEntry.split(':')
        except:
            thisSolType, pol, ant = solEntry.split(':')
            pol1 = pol
            pol2 = pol
        dir = 'pointing'

    # For RotationMeasure assuming [RotationMeasure:ant]
    elif thisSolType == 'RotationMeasure':
        dir = 'pointing'
        try:
            thisSolType, ant = solEntry.split(':')
        except:
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

    if pol1 is not None and pol2 is not None:
        if pol1 == '0' and pol2 == '0': pol = 'XX'
        if pol1 == '1' and pol2 == '0': pol = 'YX'
        if pol1 == '0' and pol2 == '1': pol = 'XY'
        if pol1 == '1' and pol2 == '1': pol = 'YY'

    if pol is not None:
        pol = re.escape(pol)
    if dir is not None:
        dir = re.escape(dir)
    if ant is not None:
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
    if parm is not None:
        parm = parm.lower()

    for st in solTabs:
        # Handle gain table separately, as we need to distinguish ampl and phase
        if solType == 'DirectionalGain' or solType == 'Gain':
            if (parm == 'ampl' or parm == 'real') and st.getType() == 'amplitude':
                if hasattr(st.obj._v_attrs, 'parmdb_type'):
                    if st.obj._v_attrs['parmdb_type'] is not None:
                        if solType in st.obj._v_attrs['parmdb_type'].split(', '):
                            solTabList.append(st)
                    else:
                        solTabList.append(st)
                else:
                    solTabList.append(st)
            elif (parm == 'phase' or parm == 'imag') and st.getType() == 'phase':
                if hasattr(st.obj._v_attrs, 'parmdb_type'):
                    if st.obj._v_attrs['parmdb_type'] is not None:
                        if solType in st.obj._v_attrs['parmdb_type'].split(', '):
                            solTabList.append(st)
                    else:
                        solTabList.append(st)
                else:
                    solTabList.append(st)
        else:
            if hasattr(st.obj._v_attrs, 'parmdb_type'):
                if st.obj._v_attrs['parmdb_type'] is not None:
                    if solType in st.obj._v_attrs['parmdb_type'].split(', '):
                        solTabList.append(st)
                else:
                    if (solType == 'RotationAngle' or solType == 'CommonRotationAngle') and st.getType() == 'rotation':
                        solTabList.append(st)
                    elif (solType == 'CommonScalarPhase' or solType == 'ScalarPhase') and st.getType() == 'scalarphase':
                        solTabList.append(st)
                    elif (solType == 'CommonScalarAmplitude' or solType == 'ScalarAmplitude') and st.getType() == 'scalaramplitude':
                        solTabList.append(st)
                    elif solType == 'Clock' and st.getType() == 'clock':
                        solTabList.append(st)
                    elif solType == 'TEC' and st.getType() == 'tec':
                        solTabList.append(st)
                    elif solType == 'RotationMeasure' and st.getType() == 'rotationmeasure':
                        solTabList.append(st)
            else:
                if (solType == 'RotationAngle' or solType == 'CommonRotationAngle') and st.getType() == 'rotation':
                    solTabList.append(st)
                elif (solType == 'CommonScalarPhase' or solType == 'ScalarPhase') and st.getType() == 'scalarphase':
                    solTabList.append(st)
                elif (solType == 'CommonScalarAmplitude' or solType == 'ScalarAmplitude') and st.getType() == 'scalaramplitude':
                    solTabList.append(st)
                elif solType == 'Clock' and st.getType() == 'clock':
                    solTabList.append(st)
                elif solType == 'TEC' and st.getType() == 'tec':
                    solTabList.append(st)
                elif solType == 'RotationMeasure' and st.getType() == 'rotationmeasure':
                    solTabList.append(st)

    if len(solTabList) == 0:
        return None
    else:
        return solTabList


def makeTECparmdb(H, solset, TECsolTab, timewidths, freq, freqwidth):
    """Returns TEC screen parmdb parameters

    H - H5parm object
    solset - solution set with TEC screen parameters
    TECsolTab = solution table with tecscreen values
    timewidths - time widths of output parmdb
    freq - frequency of output parmdb
    freqwidth - frequency width of output parmdb
    """
    global ipbar, pbar

    solset = H.getSolset(solset)

    station_dict = solset.getAnt()
    station_names = list(station_dict.keys())
    station_positions = list(station_dict.values())

    source_dict = solset.getSou()
    source_names = list(source_dict.keys())
    source_positions = list(source_dict.values())

    tec_sf = solset.getSoltab(TECsolTab)
    tec_screen, axis_vals = tec_sf.getValues()
    times = axis_vals['time']
    beta = TECsolTab._v_attrs['beta']
    r_0 = TECsolTab._v_attrs['r_0']
    height = TECsolTab._v_attrs['height']
    order = TECsolTab._v_attrs['order']
    pp = tec_sf.t.piercepoint

    N_sources = len(source_names)
    N_times = len(times)
    N_freqs = 1
    N_stations = len(station_names)
    N_piercepoints = N_sources * N_stations

    freqs = freq
    freqwidths = freqwidth
    parms = {}
    v = {}
    v['times'] = times
    v['timewidths'] = timewidths
    v['freqs'] = freqs
    v['freqwidths'] = freqwidths

    for station_name in station_names:
        for source_name in source_names:

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

            v['values'] = np.zeros((N_times, N_freqs), dtype=np.double)
            parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
            parms[parmname] = v.copy()

    for k in range(N_times):
        D = np.resize(pp[k, :, :], (N_piercepoints, N_piercepoints, 3))
        D = np.transpose(D, ( 1, 0, 2 )) - D
        D2 = np.sum(D**2, axis=2)
        C = -(D2 / (r_0**2))**(beta / 2.0) / 2.0
        tec_fit_white = np.dot(np.linalg.inv(C),
            tec_screen[:, k, :].reshape(N_piercepoints))
        pp_idx = 0
        for src, source_name in enumerate(source_names):
            for sta, station_name in enumerate(station_names):

                parmname = 'Piercepoint:X:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 0]

                parmname = 'Piercepoint:Y:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 1]

                parmname = 'Piercepoint:Z:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = pp[k, pp_idx, 2]

                parmname = 'TECfit_white:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                parmname = 'TECfit_white:0:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                parmname = 'TECfit_white:1:%s:%s' % (station_name, source_name)
                parms[parmname]['values'][k, 0] = tec_fit_white[pp_idx]

                pp_idx += 1
        pbar.update(ipbar)
        ipbar += 1

    time_start = times[0] - timewidths[0]/2
    time_end = times[-1] + timewidths[-1]/2

    v['times'] = np.array([(time_start + time_end) / 2])
    v['timewidths'] = np.array([time_end - time_start])

    v_r0 = v.copy()
    v_r0['values'] = np.array(r_0, dtype=np.double, ndmin=2)
    parms['r_0'] = v_r0

    v_beta = v.copy()
    v_beta['values'] = np.array(beta, dtype=np.double, ndmin=2)
    parms['beta'] = v_beta

    v_height = v.copy()
    v_height['values'] = np.array(height, dtype=np.double, ndmin=2)
    parms['height'] = v_height

    return parms


if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog <H5parm filename> <output globaldb/SB filename>\n'+
        _author, version='%prog '+_version.__version__)
    opt.add_option('-V', '--verbose', help='Go VeRbOsE!',
        action='store_true', default=False)
    opt.add_option('-s', '--solset', help='Name of solution set to export '
        '(default=sol000)', type='string', default='sol000')
    opt.add_option('-o', '--outfile', help='Filename of globaldb/SB to export parmdb to '
        '(default=input globaldb/SB filename)', type='string', default=None)
    opt.add_option('-r', '--root', help='Root string to prepend to input parmdb '
        'instrument directories to make the output parmdb directories '
        '(default=solution-set name)', type='string', default=None)
    opt.add_option('-t', '--soltab', help='Solution tables to export; e.g., '
        '"amplitude000, phase000" (default=all)', type='string', default='all')
    opt.add_option('-i', '--instrument', help='Name of the instrument table '
        '(default=instrument*)', type='string', default='instrument*')
    opt.add_option('-c', '--clobber', help='Clobber exising files '
        '(default=False)', action='store_true', default=False)
    (options, args) = opt.parse_args()
    global ipbar, pbar

    logger = _logging.Logger('info')
    logging = _logging.logger

    # Check options
    if len(args) != 2:
        opt.print_help()
        sys.exit()
    if options.verbose: logger.set_level("debug")

    # Check input H5parm file
    h5parmFile = args[0]
    if not os.path.exists(h5parmFile):
        logging.critical('Input H5parm file not found.')
        sys.exit(1)
    logging.info("Input H5parm filename = "+h5parmFile)

    # Open the h5parm file and get solution set names
    h5parm_in = h5parm(h5parmFile, readonly = True)
    solsetNames = h5parm_in.getSolsetNames()

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
        glob.glob(os.path.join(globaldbFile, options.instrument)) \
        if os.path.isdir(instrumentdbFile) ]
    if len(instrumentdbFiles) == 0:
        logging.critical('No parmdb table(s) found in input globaldb/SB file.')
        sys.exit(1)
    instrumentdbFiles.sort()

    # Find solution table types using the first instrumentdb
    # TODO: is there a better solution which check all the instrumentdbs?
    pdb = lofar.parmdb.parmdb(instrumentdbFiles[0])
    solTypes = list(set(x[0] for x in  (x.split(":") for x in pdb.getNames())))
    solTabs = solset.getSoltabs()
    if options.soltab != 'all':
        soltabs_to_use = [s.strip() for s in options.soltab.split(',')]
        logging.info('Using solution tables: {0}'.format(soltabs_to_use))
        solTabs_filt = []
        for s in solTabs:
            if s.name in soltabs_to_use:
                solTabs_filt.append(s)
        for s in soltabs_to_use:
            if s.name not in solset.getSoltabNames():
                logging.warning('Solution table {0} not found in input H5parm file.'.format(s.name))
        solTabs = solTabs_filt
    if len(solTabs) == 0:
        logging.critical('No solution tables found in input H5parm file')
        sys.exit(1)
    pdbSolTypes = solTypes[:]
    for solType in pdbSolTypes:
        solTabList = getSoltabFromSolType(solType, solTabs, 'ampl')
        if solTabList is None:
            # Search for type phase solutions if no ampl ones where found
            solTabList = getSoltabFromSolType(solType, solTabs, 'phase')
        if solTabList is None:
            logging.warning("Solution type {0} not found in solution set {1}. Skipping.".format(solType, solsetName))
            solTypes.remove(solType)

    # Look for tecscreen solution table in the solset. If
    # found, add to solTypes
    st_tec = None
    for st in solTabs:
        if st.getType() == 'tecscreen':
            st_tec = st
    if st_tec is not None:
        solTypes.append('TECScreen')
    solTypes = list(set(solTypes))
    logging.info('Found solution types in input parmdb and H5parm: '+', '.join(solTypes))

    # For each solType, select appropriate solutions and construct
    # the dictionary to pass to pdb.addValues()
    len_sol = {}
    for solType in solTypes:
        if solType != 'TECScreen':
            len_sol[solType] = len(pdb.getNames(solType+':*'))
        else:
            N_times = st_tec.getAxisLen(axis='time')
            len_sol[solType] = N_times

    cachedSolTabs = {}
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

        # Add default values and steps
        DefValues = pdb_in.getDefValues()
        for k, v in DefValues.items():
            pdb_out.addDefValues({k: pdb.makeDefValue(v.item(0))})
        pdb_out.setDefaultSteps(pdb_in.getDefaultSteps())

        for solType in solTypes:
            if len_sol[solType] == 0: continue

            if solType != 'TECScreen':
                solEntries = pdb_in.getNames(solType+':*')
                data = pdb_in.getValuesGrid(solType+':*')
                data_out = data.copy()
                for solEntry in solEntries:

                    pol, dir, ant, parm = parmdbToAxes(solEntry)
                    solTabList = getSoltabFromSolType(solType, solTabs, parm=parm)
                    if solTabList is None:
                        continue
                    if len(solTabList) > 1:
                        logging.warning('More than one solution table found in H5parm '
                                'matching parmdb entry "'+solType+'". Taking the first match: '+str(solTabList[0])+'.')
                    solTab = solTabList[0]

                    # search in the cache for open soltab
                    if not solTab.getType() in cachedSolTabs:
                        cachedSolTabs[solTab.getType()] = solTab
                    else:
                        solTab = cachedSolTabs[solTab.getType()]
                        solTab.setSelection()

                    freqs = data[solEntry]['freqs']
                    times = data[solEntry]['times']
                    parms = {}

                    if 'ant' in solTab.getAxesNames():
                        parms['ant'] = [ant]
                        # skip missing antennas (e.g. internationals sometimes are retained in the parmdb)
                        if not ant in solTab.getAxisValues('ant'): continue
                    if 'pol' in solTab.getAxesNames(): parms['pol'] = [pol]
                    if 'dir' in solTab.getAxesNames(): parms['dir'] = [dir]
                    if 'freq' in solTab.getAxesNames(): parms['freq'] = freqs.tolist()
                    # workaround for bbs and ndppp dealing differently with the last time slot when #timeslots%ntime != 0
                    # NDPPP has all intervals the same
                    # BBS has a maller interval in the last timeslot which is compensated here
                    if times[-1] - times[-2] < times[-2] - times[-3]: times[-1] = times[-2] + (times[-2] - times[-3])
                    if 'time' in solTab.getAxesNames(): parms['time'] = {'min':np.min(times-0.1), 'max':np.max(times+0.1)}
                    solTab.setSelection(**parms)

                    # If needed, convert Amp and Phase to Real and Imag
                    if parm == 'Real':
                        solTabList = getSoltabFromSolType(solType, solTabs, parm='phase')
                        if not solTabList[0].getType() in cachedSolTabs:
                            solTab_ph = solTabList[0]
                            cachedSolTabs[solTabList[0].getType()] = solTab_ph
                        else:
                            solTab_ph = cachedSolTabs[solTabList[0].getType()]
                        solTab_ph.setSelection(ant=[ant], pol=[pol], dir=[dir], freq=freqs.tolist(),
                            time={'min':np.min(times-0.1), 'max':np.max(times+0.1)})
                        val_phase = solTab_ph.getValues()[0]
                        val_amp = solTab.getValues()[0]
                        val = val_amp * np.cos(val_phase)
                    elif parm == 'Imag':
                        solTabList = getSoltabFromSolType(solType, solTabs, parm='ampl')
                        if not solTabList[0].getType() in cachedSolTabs:
                            solTab_amp = solTabList[0]
                            cachedSolTabs[solTabList[0].getType()] = solTab_amp
                        else:
                            solTab_amp = cachedSolTabs[solTabList[0].getType()]
                        solTab_amp.setSelection(ant=[ant], pol=[pol], dir=[dir], freq=freqs.tolist(),
                            time={'min':np.min(times-0.1), 'max':np.max(times+0.1)})
                        val_phase = solTab.getValues()[0]
                        val_amp = solTab_amp.getValues()[0]
                        val = val_amp * np.sin(val_phase)
                    else:
                        val = solTab.getValues()[0]

                    # Apply flags
                    weights = solTab.getValues(weight=True)[0]

                    # etienne part; if it is borken, curse his name
                    # check whether this is clock or tec; if so, reshape properly to account for all freqs in the parmdb
                    # anyway these tables are freq-indep
                    if solType == "Clock":# or solType == "TEC" or solType == "RotationMeasure":
                        # find freq-dimensionality 
                        nfreq = freqs.shape[0]
                        #print val.shape
                        # reshape such that all freq arrays are filled properly
                        val = np.tile( val, np.append([nfreq], np.ones(len(val.shape),dtype=np.int) ) )
                        #print val.shape
                        weights = np.tile( weights, np.append([nfreq], np.ones(len(weights.shape),dtype=np.int) ) )

                    flags = np.zeros(shape=weights.shape, dtype=bool)
                    flags[np.where(weights == 0)] = True
                    if parm == 'Real':
                        weights2 = solTab_ph.getValues(weight=True)[0]
                        flags[np.where(weights2 == 0)] = True
                    if parm == 'Imag':
                        weights2 = solTab_amp.getValues(weight=True)[0]
                        flags[np.where(weights2 == 0)] = True
                    np.putmask(val, flags, np.nan)

                    shape = data_out[solEntry]['values'].shape
                    #print "shape"
                    #print 'parmdb', shape
                    #print 'h5parm', val.shape
                    #print sf.getAxesNames()
                    #print sf.getType()
                    #print "parmdb", times
                    #for t in times: print '%.1f' % t
                    #print "h5parm", sf.time
                    #for t in sf.time: print '%.1f' % t
                    try:
                        data_out[solEntry]['values'] = val.T.reshape(shape)
                    except ValueError as err:
                        logging.critical('Mismatch between parmdb table and H5parm '
                        'solution table: Differing number of frequencies and/or times')
                        sys.exit(1)
                ipbar += 1
                pbar.update(ipbar)
            else:
                # Handle TECScreen parmdb
                #
                # Get timewidths, freqwidth and freq from first (non-TEC, phase)
                # solentry
                for nonTECsolType in pdbSolTypes:
                    if nonTECsolType != 'TECScreen' and 'Phase' in nonTECsolType:
                        break
                parmname = pdb_in.getNames(nonTECsolType+':*')[0]
                timewidths = pdb_in.getValuesGrid(parmname)[parmname]['timewidths']
                freqwidth = pdb.getValuesGrid(parmname)[parmname]['freqwidths'][0]
                freq = pdb.getValuesGrid(parmname)[parmname]['freqs'][0]
                data_out = makeTECparmdb(h5parm_in, solset, st_tec, timewidths, freq, freqwidth)

            pdb_out.addValues(data_out)

        pbar.finish()

    logging.info('Done.')
