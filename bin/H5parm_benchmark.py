#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool comapre the performances of H5parm with parmdb.

# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)"

import sys, os, time
import numpy as np
import logging
import lofar.parmdb
from losoto import _version
from losoto import _logging
from losoto.h5parm import h5parm, solFetcher, solWriter

# Options
import optparse
opt = optparse.OptionParser(usage='%prog -p H5parm [-s solset] -g parmdb [-n 100]\n'\
                +_author, version='%prog '+_version.__version__)
opt.add_option('-p', '--h5parm', help='H5parm name', type='string', default='')
opt.add_option('-g', '--parmdb', help='Parmdb name', type='string', default='')
opt.add_option('-s', '--solset', help='Solution-set name (default=sol000)', type='string', default='sol000')
opt.add_option('-n', '--numiter', help='Number of iterations (default=100)', type=int, default=100)
(options, args) = opt.parse_args()

_logging.setLevel('debug')

n = options.numiter

solset = options.solset
h5parmFile = options.h5parm
H5 = h5parm(h5parmFile, readonly=False)
H = solFetcher(H5.getSoltab(solset,'rotation000'))
logging.info("H5parm filename = "+h5parmFile)

parmdbFile = options.parmdb
P = lofar.parmdb.parmdb(parmdbFile)
P2 = lofar.parmdb.parmdb('tmp.parmdb', create=True)
logging.info("parmdb filename = "+parmdbFile)

######################################################
logging.info("### Read all frequencies for a pol/dir/station")

start = time.clock()
for i in xrange(n):
    Pfreqs = P.getValuesGrid('RotationAngle:CS001LBA:3C196')['RotationAngle:CS001LBA:3C196']['freqs']
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
for i in xrange(n):
    H.setSelection(dir='3C196',ant='CS001LBA')
    Hfreqs = H.freq
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

#print "Equal?",(Pfreqs == Hfreqs).all()

######################################################
logging.info("### Read all times for a pol/dir/station")

start = time.clock()
for i in xrange(n):
    Ptimes = P.getValuesGrid('RotationAngle:CS001LBA:3C196')['RotationAngle:CS001LBA:3C196']['times']
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
for i in xrange(n):
    H.setSelection(dir='3C196',ant='CS001LBA')
    Htimes = H.time
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

#print "Equal?",(Ptimes == Htimes).all()

######################################################
logging.info("### Read all rotations for 1 station (slice in time)")

start = time.clock()
for i in xrange(n):
    Prot = P.getValuesGrid('CommonRotationAngle:CS001LBA',stime=Ptimes[30],etime=Ptimes[-30])['CommonRotationAngle:CS001LBA']['values'].transpose()[0]
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
for i in xrange(n):
    H.setSelection(dir='pointing', ant='CS001LBA', time={'min':Htimes[31],'max':Htimes[-29]})
    Hrot = H.getValues(retAxesVals = False)
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

#print "Equal?", (Prot == np.squeeze(Hrot)).all()

######################################################
logging.info("### Read all rotations for all station (slice in time)")

start = time.clock()
for i in xrange(n):
    Protd = P.getValuesGrid('CommonRotationAngle:*',stime=Ptimes[30],etime=Ptimes[-30])
    # construct the matrix
    Prot = []
    for ant in sorted(Protd.iterkeys()):
        Prot.append(Protd[ant]['values'].transpose()[0])
    Prot = np.array(Prot)
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
for i in xrange(n):
    H.setSelection(dir='pointing', time={'min':Htimes[31],'max':Htimes[-29]})
    Hrot = H.getValues(retAxesVals = False)
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

#print "Equal?", (Prot == np.squeeze(Hrot)).all()

######################################################
logging.info("### Read all rotations for remote stations (slice in ant)")

start = time.clock()
for i in xrange(n):
    Protd = P.getValuesGrid('CommonRotationAngle:RS*')
    # construct the matrix
    Prot = []
    for ant in sorted(Protd.iterkeys()):
        Prot.append(Protd[ant]['values'].transpose()[0])
    Prot = np.array(Prot)
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
for i in xrange(n):
    H.setSelection(dir='pointing', ant='RS*')
    Hrot = H.getValues(retAxesVals = False)
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

#print "Equal?", (Prot == np.squeeze(Hrot)).all()

######################################################
# read+write
logging.info("### Read all rotations for a dir/station and write them back")
Hw = solWriter(H5.getSoltab(solset,'rotation000'))

start = time.clock()
for i in xrange(n):
    Prot = P.getValuesGrid('CommonRotationAngle:CS001LBA')
    Prot = {'test'+str(i):Prot['CommonRotationAngle:CS001LBA']}
    P2.addValues(Prot)
    # parmdb write?
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
for i in xrange(n):
    H.setSelection(dir='pointing', ant='CS001LBA')
    Hrot = H.getValues(retAxesVals = False)
    Hw.setSelection(dir='pointing', ant='CS001LBA')
    Hw.setValues(Hrot)
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

######################################################
# read whole file
logging.info("### Read and tabulate the whole file")

start = time.clock()
val = P.getValuesGrid('')
elapsed = (time.clock() - start)
logging.info("parmdb -- "+str(elapsed)+" s.")

start = time.clock()
H.setSelection()
val, axes = H.getValues()
#print "Shape:", val.shape
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")




del H
del P
del P2
os.system('rm -r tmp.parmdb')
