#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool uses a gds file to collect all the necessary fiels from the
# cluster and create a set of parmdbs that can be used by H5parm_creator.py
# to create an H5parm file.
# It also collect the necessary sky, FIELD and ANTENNA tables from one MS.

# Authors:
# Francesco de Gasperin
# Bas van der Tol
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)\n\
Bas van der Tol (vdtol@strw.leidenuniv.nl)"

import sys, os, re
import logging
import lofar.parameterset
import losoto._version
import losoto._logging

def splitgds(gdsFile, wd='', id='part'):
    """Split gds file in multiple files
    """
    ps = lofar.parameterset.parameterset(gdsFile)
    clusterdesc = ps.getString('ClusterDesc')
    starttime = ps.getString('StartTime')
    endtime = ps.getString('EndTime')
    steptime = ps.getString('StepTime')
    N = ps.getInt('NParts')
    gdsList = []
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
        gdsList.append(partname)
    return gdsList


if __name__=='__main__':
    # Options
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v] [-o] [-d gds] [-g output globaldb] \n'\
            +_author, version='%prog '+losoto._version.__version__)
    opt.add_option('-v', '--verbose', help='Go VeRbOsE! (default=False)', action='store_true', default=False)
    opt.add_option('-o', '--overwrite', help='Overwrite an existing globaldb (default=False)', action='store_true', default=False)
    opt.add_option('-d', '--gds', help='Gds file used to construct the globaldb', type='string', default='')
    opt.add_option('-g', '--globaldb', help='Output globaldb name (default=globaldb)', type='string', default='globaldb')
    (options, args) = opt.parse_args()

    # Check options
    if len(args) != 0:
        opt.print_help()
        sys.exit()
    if options.verbose: losoto._logging.setLevel('debug')

    overwrite = options.overwrite
    gdsFile = options.gds
    logging.info("GDS filename = "+gdsFile)
    globaldbFile = options.globaldb
    logging.info("globaldb filename = "+globaldbFile)

    # hardcoded, maybe can be turned in an option...
    instrumentName='instrument'

    if not os.path.isfile(gdsFile):
        logging.critical('gdsFile '+gdsFile+' not found.')
        sys.exit()

    if os.path.exists(globaldbFile):
        if overwrite:
            os.system('rm -rf '+globaldbFile)
            os.makedirs(globaldbFile)
        else:
            logging.warning(globaldbFile+' already exists. I will not overwrite existing files.')
    else:
        os.makedirs(globaldbFile)

    # Create an instrumentdb named as gdsFile.instrumentName
    # which is like the gds file but points to the "instrument" parmdb table
    # inside each MS.

    instrumentdbFile = os.path.join(globaldbFile, \
                os.path.splitext(os.path.basename(gdsFile))[0] \
                + os.path.extsep + instrumentName)
    p = re.compile('(^Part\\d*.FileName\\s*=\\s*\\S*)')
    gdsFileR = open(gdsFile)
    instrumentdbFileW = open(instrumentdbFile, 'w')
    instrumentdbFileW.writelines([p.sub('\\1%s%s' % (os.path.sep,\
            instrumentName), l) for l in gdsFileR.readlines()])
    gdsFileR.close()
    instrumentdbFileW.close()

    # split the gdsFile and the instrumentdbFile for each SB
    gdsFiles = splitgds(gdsFile, wd=globaldbFile, id='part')
    instrumentdbGdsFiles = splitgds(instrumentdbFile, wd=globaldbFile, id='instrument')

    # Collect all the instrument tables
    instrumentdbFiles = []
    for instrumentdbGdsFile in instrumentdbGdsFiles:
        instrumentdbParset = lofar.parameterset.parameterset(instrumentdbGdsFile)
        instrumentdbRemoteFile = instrumentdbParset.getString('Part0.FileName')
        instrumentdbHostname = instrumentdbParset.getString('Part0.FileSys').split(':')[0]
        instrumentdbFile = os.path.splitext(instrumentdbGdsFile)[0]
        if not os.path.exists(instrumentdbFile):
            logging.info("Collecting "+instrumentdbFile)
            if instrumentdbHostname == 'localhost':
                os.system('cp -r %s %s > /dev/null' % (instrumentdbRemoteFile, instrumentdbFile))
            else:
                os.system('scp -r %s:%s %s > /dev/null' % (instrumentdbHostname, instrumentdbRemoteFile, instrumentdbFile))
        else:
            logging.info("Skipping "+instrumentdbFile)
        instrumentdbFiles.append(instrumentdbFile)

    gdsParset = lofar.parameterset.parameterset(gdsFiles[0])
    # instrumentdbParset =

    hostname = gdsParset.getString('Part0.FileSys').split(':')[0]
    msname = gdsParset.getString('Part0.FileName')
    # Collect the skydb from the first SB
    skydbFile = os.path.join(globaldbFile, 'sky')
    if not os.path.exists(skydbFile):
        logging.info("Collecting the skydb")
        if hostname == 'localhost':
            os.system('cp -r %s/sky %s > /dev/null' % (msname, skydbFile))
        else:
            os.system('scp -r %s:%s/sky %s > /dev/null' % (hostname, msname, skydbFile))
    # Collect the ANTENNA table from the first SB
    antennaFile = os.path.join(globaldbFile, 'ANTENNA')
    if not os.path.exists(antennaFile):
        logging.info("Collecting the antenna table")
        if hostname == 'localhost':
            os.system('cp -r %s/ANTENNA %s > /dev/null' % (msname, antennaFile))
        else:
            os.system('scp -r %s:%s/ANTENNA %s > /dev/null' % (hostname, msname, antennaFile))
    # Collect the FILED table from the first SB
    fieldFile = os.path.join(globaldbFile, 'FIELD')
    if not os.path.exists(fieldFile):
        logging.info("Collecting the field table")
        if hostname == 'localhost':
            os.system('cp -r %s/FIELD %s > /dev/null' % (msname, fieldFile))
        else:
            os.system('scp -r %s:%s/FIELD %s > /dev/null' % (hostname, msname, fieldFile))

    logging.info("Done.")
