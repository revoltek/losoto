#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, time, re
from math import *

import numpy
import scipy.optimize

import lofar.parmdb
import lofar.parameterset
import pyrap.tables as pt
import tables


def make_instrumentdb(gdsfilename, instrument_name, globaldb):
    """ Create an instrumentdb named as gdsfilename.instrument_name
    which is like the gds file but points to the "instrument" parmdb tables
    inside each MS.
    """
    instrumentdb_name = os.path.join(globaldb,
            os.path.splitext(os.path.basename(gdsfilename))[0]
            + os.path.extsep + instrument_name)
    p = re.compile('(^Part\\d*.FileName\\s*=\\s*\\S*)')
    gdsfile = open(gdsfilename)
    instrumentdb_file = open(instrumentdb_name, 'w')
    instrumentdb_file.writelines([p.sub('\\1%s%s' % (os.path.sep,
                                 instrument_name), l) for l in
                                 gdsfile.readlines()])
    gdsfile.close()
    instrumentdb_file.close()
    return instrumentdb_name


def splitgds(gdsname, wd='', id='part'):
    """Split gds file in multiple files
    """
    ps = lofar.parameterset.parameterset(gdsname)
    clusterdesc = ps.getString('ClusterDesc')
    starttime = ps.getString('StartTime')
    endtime = ps.getString('EndTime')
    steptime = ps.getString('StepTime')
    N = ps.getInt('NParts')
    gdslist = []
    for i in range(N):
        partname = os.path.join(wd, '%s-%i.gds' % (id, i))
        ps_part = ps.makeSubset('Part%i.' % i, 'Part0.')
        NChan = ps_part.getString('Part0.NChan')
        StartFreqs = ps_part.getString('Part0.StartFreqs')
        EndFreqs = ps_part.getString('Part0.EndFreqs')

        ps_part.add('Name', os.path.basename(partname))
        ps_part.add('ClusterDesc', clusterdesc)
        ps_part.add('StartTime', starttime)
        ps_part.add('EndTime', endtime)
        ps_part.add('StepTime', steptime)
        ps_part.add('NChan', NChan)
        ps_part.add('StartFreqs', StartFreqs)
        ps_part.add('EndFreqs', EndFreqs)
        ps_part.add('NParts', '1')
        ps_part.writeFile(partname)
        gdslist.append(partname)
    return gdslist


def get_station_list(pdb, station_pattern_list, DirectionalGainEnable):
    station_list = []
    for pattern in station_pattern_list:
        parmname_list = pdb.getNames({True: 'DirectionalGain:?:?:*:'
                + pattern + ':*', False: 'Gain:?:?:*:'
                + pattern}[DirectionalGainEnable])
        station_list.extend(sorted(set([n.split(':')[{True: -2,
                            False: -1}[DirectionalGainEnable]] for n in
                            parmname_list])))
    return station_list


def fillarray(a, v):
    print a.shape, a.chunkshape
    for idx in product(*[xrange(0, s, c) for (s, c) in zip(a.shape,
                       a.chunkshape)]):
        s = tuple([slice(i, min(i + c, s)) for (i, s, c) in zip(idx,
                  a.shape, a.chunkshape)])
        a[s] = v


def product(*args, **kwds):
    """
    product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    """
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x + [y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)


class ClassMakeHDF5:

    def __init__(self):
        self.clusterdesc = ''
        self.sky_name = 'sky'
        self.instrument_name = 'instrument'
        self.sources = []
        self.stations = []
        self.globaldb = 'globaldb'
        self.PhasorsEnable = False
        self.GainEnable = True  # False
        self.polarizations = [0, 1]
        self.N_pol = len(self.polarizations)
        self.DirectionalGainEnable = False

    def load_globaldb(self, globaldb):
        self.globaldb = globaldb
        self.hdf5 = tables.openFile(os.path.join(globaldb,'ionmodel.hdf5'), 'r+')

        self.stations = self.hdf5.root.stations.cols.name
        self.station_positions = self.hdf5.root.stations.cols.position
        self.array_center = self.hdf5.root.array_center
        self.N_stations = len(self.stations)

        self.sources = self.hdf5.root.sources[:]['name']
        self.source_positions = self.hdf5.root.sources[:]['position']
        self.N_sources = len(self.sources)
        self.N_piercepoints = self.N_sources * self.N_stations

        self.pointing = self.hdf5.root.pointing

        self.freqs = self.hdf5.root.freqs
        self.polarizations = self.hdf5.root.polarizations
        self.N_pol = len(self.polarizations)

        self.flags = self.hdf5.root.flags

        for varname in [
            'amplitudes',
            'phases',
            'Clock',
            'TEC',
            'TECfit',
            'TECfit_white',
            'offsets',
            'times',
            'timewidths',
            'piercepoints',
            'facets',
            'facet_piercepoints',
            'n_list',
            'STEC_facets',
            ]:
            if varname in self.hdf5.root:
                self.__dict__.update([(varname,
                        self.hdf5.getNode(self.hdf5.root, varname))])

        self.N_stations = len(self.stations)
        self.N_sources = len(self.sources)

    def load_gds(
        self,
        # TODO: why is this a list?
        gdsfiles,
        clusterdesc,
        globaldb='globaldb',
        sky_name='sky',
        instrument_name='instrument',
        stations=[],
        sources=[],
        ):

        self.gdsfiles = gdsfiles
        self.instrument_name = instrument_name
        self.globaldb = globaldb

        if not os.path.exists(globaldb):
            os.makedirs(globaldb)

        self.instrumentdb_name_list = []
        for gdsfile in gdsfiles:
            instrumentdb_name = os.path.splitext(gdsfile)[0] \
                + os.path.extsep + instrument_name
            if not os.path.exists(instrumentdb_name):
                instrumentdb_name = make_instrumentdb(gdsfile,
                        instrument_name, globaldb)
            self.instrumentdb_name_list.append(instrumentdb_name)

        gdsfiles = []
        for (idx, gdsfile) in zip(range(len(self.gdsfiles)), self.gdsfiles):
            gdsfiles.extend(splitgds(gdsfile, wd=self.globaldb, id='part-%i' % idx))
        self.gdsfiles = gdsfiles

        instrumentdb_name_list = []
        for (idx, instrumentdb_name) in \
            zip(range(len(self.instrumentdb_name_list)),
                self.instrumentdb_name_list):
            instrumentdb_name_list.extend(splitgds(instrumentdb_name,
                    wd=self.globaldb, id='instrument-%i' % idx))
        self.instrumentdb_name_list = instrumentdb_name_list

        instrumentdb_name_list = []
        for instrumentdb_name in self.instrumentdb_name_list:
            instrumentdb_parset = lofar.parameterset.parameterset(instrumentdb_name)
            instrumentdbfilename = instrumentdb_parset.getString('Part0.FileName')
            instrumentdbhostname = instrumentdb_parset.getString('Part0.FileSys').split(':')[0]
            instrumentdb_name = os.path.splitext(instrumentdb_name)[0]
            if not os.path.exists(instrumentdb_name):
                os.system('scp -r %s:%s %s' % (instrumentdbhostname,
                          instrumentdbfilename, instrumentdb_name))
            instrumentdb_name_list.append(instrumentdb_name)
        self.instrumentdb_name_list = instrumentdb_name_list

        gdsfile = gdsfiles[0]
        instrumentdb = lofar.parmdb.parmdb(self.instrumentdb_name_list[0])

        self.hdf5 = tables.openFile(os.path.join(globaldb, 'ionmodel.hdf5'), 'w')

        gdsparset = lofar.parameterset.parameterset(gdsfile)

        skydbfilename = os.path.join(gdsparset.getString('Part0.FileName'), sky_name)
        skydbhostname = gdsparset.getString('Part0.FileSys').split(':')[0]
        skydbname = globaldb + '/sky'
        if not os.path.exists(skydbname):
            os.system('scp -r %s:%s %s' % (skydbhostname,
                      skydbfilename, skydbname))
        skydb = lofar.parmdb.parmdb(skydbname)

        gdsparset = lofar.parameterset.parameterset(gdsfile)
        msname = gdsparset.getString('Part0.FileName')
        mshostname = gdsparset.getString('Part0.FileSys').split(':')[0]
        antenna_table_name = os.path.join(globaldb, 'ANTENNA')
        if not os.path.exists(antenna_table_name):
            os.system('scp -r %s:%s/ANTENNA %s' % (mshostname, msname, antenna_table_name))
        field_table_name = os.path.join(globaldb, 'FIELD')
        if not os.path.exists(field_table_name):
            os.system('scp -r %s:%s/FIELD %s' % (mshostname, msname, field_table_name))

        if len(stations) == 0:
            stations = ['*']
        self.stations = get_station_list(instrumentdb, stations, self.DirectionalGainEnable)
        self.N_stations = len(self.stations)

        antenna_table = pt.table(globaldb + '/ANTENNA')
        name_col = antenna_table.getcol('NAME')
        position_col = antenna_table.getcol('POSITION')
        self.station_positions = [position_col[name_col.index(station_name)]
             for station_name in self.stations]
        antenna_table.close()

        station_table = self.hdf5.createTable(self.hdf5.root, 'stations'
                , {'name': tables.StringCol(40),
                'position': tables.Float64Col(3)})
        row = station_table.row
        for (station, position) in zip(self.stations, self.station_positions):
            row['name'] = station
            row['position'] = position
            row.append()
        station_table.flush()

        self.array_center = numpy.array(self.station_positions).mean(axis=0).tolist()
        self.hdf5.createArray(self.hdf5.root, 'array_center', self.array_center)

        field_table = pt.table(globaldb + '/FIELD')
        phase_dir_col = field_table.getcol('PHASE_DIR')
        self.pointing = phase_dir_col[0, 0, :]
        field_table.close()
        self.hdf5.createArray(self.hdf5.root, 'pointing', self.pointing)

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

        source_table = self.hdf5.createTable(self.hdf5.root, 'sources',
                {'name': tables.StringCol(40),
                'position': tables.Float64Col(2)})
        row = source_table.row
        for (source, position) in zip(self.sources, self.source_positions):
            row['name'] = source
            row['position'] = position
            row.append()
        source_table.flush()

        if self.PhasorsEnable:
            infix = ('Ampl', 'Phase')
        else:
            infix = ('Real', 'Imag')

        if self.GainEnable:
            parmname0 = ':'.join(['Gain', str(self.polarizations[0]),
                                 str(self.polarizations[0]), infix[1],
                                 self.stations[0]])
            v0 = instrumentdb.getValuesGrid(parmname0)[parmname0]
        if self.DirectionalGainEnable:
            parmname0 = ':'.join([
                'DirectionalGain',
                str(self.polarizations[0]),
                str(self.polarizations[0]),
                infix[1],
                self.stations[0],
                self.sources[0],
                ])
            v0 = instrumentdb.getValuesGrid(parmname0)[parmname0]

        self.freqs = []
        self.freqwidths = []

        # First collect all frequencies
        # We need them beforehand to sort the frequencies (frequencies are not necessarily given in sorted order)

        for instrumentdb_name in self.instrumentdb_name_list:
            instrumentdb = lofar.parmdb.parmdb(instrumentdb_name)
            v0 = instrumentdb.getValuesGrid(parmname0)[parmname0]
            freqs = v0['freqs']
            self.freqs = numpy.concatenate([self.freqs, freqs])
            self.freqwidths = numpy.concatenate([self.freqwidths, v0['freqwidths']])

        # Sort frequencies, find both the forward and inverse mapping
        # Mappings are such that
        #    sorted_freqs = unsorted_freqs[sorted_freq_idx]
        #    sorted_freqs[inverse_sorted_freq_idx] = unsorted_freqs
        # We will use the following form
        #    sorted_freqs[inverse_sorted_freq_idx[selection]] = unsorted_freqs[selection]
        # to process chunks (=selections) of unsorted data and store them in sorted order

        sorted_freq_idx = sorted(range(len(self.freqs)),
                                 key=lambda idx: self.freqs[idx])
        inverse_sorted_freq_idx = sorted(range(len(self.freqs)),
                key=lambda idx: sorted_freq_idx[idx])

        self.freqs = self.freqs[sorted_freq_idx]
        self.freqwidths = self.freqwidths[sorted_freq_idx]
        self.hdf5.createArray(self.hdf5.root, 'freqs', self.freqs)
        self.N_freqs = len(self.freqs)

        self.times = v0['times']
        self.timewidths = v0['timewidths']
        self.hdf5.createArray(self.hdf5.root, 'times', self.times)
        self.hdf5.createArray(self.hdf5.root, 'timewidths', self.timewidths)
        self.N_times = len(self.times)

        self.hdf5.createArray(self.hdf5.root, 'polarizations', self.polarizations)

        chunkshape = ( 1024, 32, 1, 1, 1, 1, )
        self.phases = self.hdf5.createCArray(self.hdf5.root, 'phases',
                tables.Float32Atom(), shape=(  # shape=(self.N_times, self.N_freqs, self.N_stations, self.N_sources, self.N_pol),
                self.N_times,
                self.N_freqs,
                self.N_stations,
                self.N_sources,
                2, 2 ), chunkshape=chunkshape)
        self.amplitudes = self.hdf5.createCArray(self.hdf5.root,
                'amplitudes', tables.Float32Atom(), shape=(  # shape=(self.N_times, self.N_freqs, self.N_stations, self.N_sources, self.N_pol),
                 self.N_times,
                 self.N_freqs,
                 self.N_stations,
                 self.N_sources,
                 2, 2 ), chunkshape=chunkshape)
        fillarray(self.amplitudes, 1.0)
        self.flags = self.hdf5.createCArray(self.hdf5.root, 'flags',
                tables.Float32Atom(), shape=(self.N_times,self.N_freqs))

        freq_idx = 0
        for (gdsfile, instrumentdb_name, gdsfile_idx) in zip(gdsfiles,
                self.instrumentdb_name_list, range(len(gdsfiles))):
            print 'Reading %s (%i/%i)' % (gdsfile, gdsfile_idx + 1,
                    len(gdsfiles))

            instrumentdb = lofar.parmdb.parmdb(instrumentdb_name)
            v0 = instrumentdb.getValuesGrid(parmname0)[parmname0]
            freqs = v0['freqs']
            N_freqs = len(freqs)
            sorted_freq_selection = inverse_sorted_freq_idx[freq_idx:freq_idx + N_freqs]

            try:
                self.flags[:, sorted_freq_selection] = \
                    instrumentdb.getValuesGrid('flags')['flags']['values']
            except KeyError:
                pass

            # print "zip:",zip(self.polarizations, range(len(self.polarizations)))
            # for pol, pol_idx in zip(self.polarizations, range(len(self.polarizations))):

            for pols in [(i, j) for i in range(2) for j in range(2)]:
                (pol0, pol1) = pols
                print 'pol: %s, pol: %s ' % (str(pol0), str(pol1))
                for (station, station_idx) in zip(self.stations,
                        range(len(self.stations))):
                    if self.GainEnable:
                        parmname0 = ':'.join(['Gain', str(pol0),
                                str(pol1), infix[0], station])
                    parmname1 = ':'.join(['Gain', str(pol0), str(pol1),
                            infix[1], station])
                    print 'parameter:', parmname0, parmname1
                    if self.PhasorsEnable:
                        gain_phase = instrumentdb.getValuesGrid(parmname1)[parmname1]['values']
                        self.phases[:, sorted_freq_selection,
                                    station_idx, :, pol0, pol1] = \
                            resize(gain_phase, (self.N_sources,
                                   N_freqs, self.N_times)).T
                        try:
                            gain_amplitude = instrumentdb.getValuesGrid(parmname0)[parmname0]['values']
                        except KeyError:
                            self.amplitudes[:, sorted_freq_selection,
                                    station_idx, :, pol0, pol1] = \
                                numpy.ones((self.N_times, N_freqs,
                                    self.N_sources))
                        else:
                            self.amplitudes[:, sorted_freq_selection,
                                    station_idx, :, pol0, pol1] = \
                                numpy.resize(gain_amplitudes,
                                    (self.N_sources, N_freqs,
                                    self.N_times)).T
                    else:
                        gain_real = \
                            instrumentdb.getValuesGrid(parmname0)[parmname0]
                        gain_imag = \
                            instrumentdb.getValuesGrid(parmname1)[parmname1]
                        self.phases[:, sorted_freq_selection,
                                    station_idx, :, pol0, pol1] = \
                            numpy.resize(numpy.arctan2(gain_imag['values'
                                ], gain_real['values']),
                                (self.N_sources, N_freqs,
                                self.N_times)).T
                        self.amplitudes[:, sorted_freq_selection,
                                station_idx, :, pol0, pol1] = \
                            numpy.resize(numpy.sqrt(gain_imag['values']
                                ** 2 + gain_real['values'] ** 2),
                                (self.N_sources, N_freqs,
                                self.N_times)).T
                if self.DirectionalGainEnable:
                    for (source, source_idx) in zip(self.sources,
                            range(len(self.sources))):
                        parmname0 = ':'.join([
                            'DirectionalGain',
                            str(pol0),
                            str(pol1),
                            infix[0],
                            station,
                            source])
                        parmname1 = ':'.join([
                            'DirectionalGain',
                            str(pol0),
                            str(pol1),
                            infix[1],
                            station,
                            source])
                        if self.PhasorsEnable:
                            gain_phase = \
                                instrumentdb.getValuesGrid(parmname1)[parmname1]['values']
                            self.phases[:, sorted_freq_selection,
                                    station_idx, source_idx, pol0,
                                    pol1] += gain_phase
                            try:
                                gain_amplitude = \
                                    instrumentdb.getValuesGrid(parmname0)[parmname0]['values']
                            except KeyError:
                                pass
                            else:
                                self.amplitudes[:,
                                        sorted_freq_selection,
                                        station_idx, source_idx, pol0,
                                        pol1] *= gain_amplitude
                        else:
                            gain_real = instrumentdb.getValuesGrid(parmname0)[parmname0]['values']
                            gain_imag = instrumentdb.getValuesGrid(parmname1)[parmname1]['values']
                            l = min(gain_real.shape[0],
                                    gain_imag.shape[0],
                                    self.phases.shape[0])
                            gain_real = gain_real[0:l, :]
                            gain_imag = gain_imag[0:l, :]
                            self.phases[0:l, sorted_freq_selection,
                                    station_idx, source_idx, pol0,
                                    pol1] += numpy.arctan2(gain_imag,
                                    gain_real)
                            self.amplitudes[0:l, sorted_freq_selection,
                                    station_idx, source_idx, pol0,
                                    pol1] *= numpy.sqrt(gain_real ** 2
                                    + gain_imag ** 2)
            freq_idx += N_freqs

        if self.flags.shape != self.phases.shape[0:2]:
            self.flags = numpy.zeros(self.phases.shape[0:2])


