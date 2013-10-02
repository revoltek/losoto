#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool comapre the performances of H5parm with parmdb.
# ./losoto_test.py -p ../examples/global-comp.h5 -g ../examples/L99289-cal_SB081.MS/instrument/ -s sol000

# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)"

import sys, os, time
import numpy as np
import logging
import lofar.parmdb
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm, solFetcher

# Options
import optparse
opt = optparse.OptionParser(usage='%prog [-v] -p H5parm [-s solset] -g parmdb \n'\
                +_author, version='%prog '+losoto._version.__version__)
opt.add_option('-p', '--h5parm', help='H5parm name', type='string', default='')
opt.add_option('-g', '--parmdb', help='Parmdb name', type='string', default='')
opt.add_option('-s', '--solset', help='Solution-set name (default=sol000)', type='string', default='sol000')
(options, args) = opt.parse_args()

losoto._logging.setVerbose()

solset = options.solset
h5parmFile = options.h5parm
H5 = h5parm(h5parmFile)
H = solFetcher(H5.getSoltab(solset,'amplitude000'))
H.makeSelection(dir='3C196',ant='CS001LBA',pol='XX')
H2 = solFetcher(H5.getSoltab(solset,'phase000'))
H2.makeSelection(dir='3C196',ant='CS001LBA',pol='XX')
logging.info("H5parm filename = "+h5parmFile)
parmdbFile = options.parmdb
P = lofar.parmdb.parmdb(parmdbFile)
logging.info("parmdb filename = "+parmdbFile)

######################################################
# read frequencies
logging.info("### Read all frequencies for a pol/dir/station")

start = time.clock()
Pfreqs = P.getValuesGrid('DirectionalGain:0:0:Real:CS001LBA:3C196')['DirectionalGain:0:0:Real:CS001LBA:3C196']['freqs']
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
Hfreqs = H.freq
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

print "Equal?",(Pfreqs == Hfreqs).all()

######################################################
# read times
logging.info("### Read all times for a pol/dir/station")

start = time.clock()
Ptimes = P.getValuesGrid('DirectionalGain:0:0:Real:CS001LBA:3C196')['DirectionalGain:0:0:Real:CS001LBA:3C196']['times']
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
H.makeSelection(dir='3C196',ant='CS001LBA',pol='XX')
Htimes = H.time
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

print Htimes
print "Equal?",(Ptimes == Htimes).all()

######################################################
## read amp solutions
#logging.info("### Read all amplitudes for a pol/dir/station")
#
#start = time.clock()
#Preal = P.getValuesGrid('DirectionalGain:0:0:Real:CS001LBA:3C196')['DirectionalGain:0:0:Real:CS001LBA:3C196']['values']
#Pim = P.getValuesGrid('DirectionalGain:0:0:Imag:CS001LBA:3C196')['DirectionalGain:0:0:Imag:CS001LBA:3C196']['values']
#Pamp = np.sqrt((Preal**2)+(Pim**2)).transpose()[0]
#elapsed = (time.clock() - start)
#logging.info("PARMDB -- "+str(elapsed)+" s.")
#
#start = time.clock()
#Hamp = H.val
#elapsed = (time.clock() - start)
#logging.info("H5parm -- "+str(elapsed)+" s.")
#
#print "Equal?", (Pamp == Hamp).all()
#
######################################################
## read real solutions
#logging.info("### Read all real for a pol/dir/station")
#
#start = time.clock()
#Preal = P.getValuesGrid('DirectionalGain:0:0:Real:CS001LBA:3C196')['DirectionalGain:0:0:Real:CS001LBA:3C196']['values'].transpose()[0]
#elapsed = (time.clock() - start)
#logging.info("PARMDB -- "+str(elapsed)+" s.")
#
#start = time.clock()
#Hamp = H.val
#Hphase = H2.val
#Hreal = (Hamp * (np.cos(Hphase) + np.sin(Hphase)*1j)).real
#elapsed = (time.clock() - start)
#logging.info("H5parm -- "+str(elapsed)+" s.")
#
#print "Equal?", (Preal == Hreal).all()

######################################################
# read a slice in time
logging.info("### Read all rotations for a dir/station (slice in time)")

start = time.clock()
Prot = P.getValues('CommonRotationAngle:CS001LBA',stime=4.868879011e+09,etime=4.86888401e+09)['CommonRotationAngle:CS001LBA']['values'].transpose()[0]
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
Hrot = [x['val'] for x in H.t.where("(time > 4.86887901e+09) & (time < 4.86888401e+09) & (ant == 'CS001LBA') & (dir == 'pointing')")]
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

print "Equal?", (Prot == Hrot).all()

######################################################
# read rotation solutions (== no conversions)
logging.info("### Read all rotations for a dir/station")

start = time.clock()
Prot = P.getValuesGrid('CommonRotationAngle:CS001LBA')['CommonRotationAngle:CS001LBA']['values'].transpose()[0]
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

H = solFetcher(H5.getSoltab(solset,'rotation000'))
H.makeSelection(ant='CS001LBA', dir='pointing')

start = time.clock()
Hrot = H.val[:,0,0,0]
elapsed = (time.clock() - start)
logging.info("H5parm H.val -- "+str(elapsed)+" s.")

print "Equal?", (Prot == Hrot).all()

start = time.clock()
Hrot2 = np.array( [x['val'] for x in H.getRowsIterator()] )
elapsed = (time.clock() - start)
logging.info("H5parm getRowsIterator() -- "+str(elapsed)+" s.")

print "Equal?", (Prot == Hrot2).all()

start = time.clock()
Hrot3, axes = H.getValuesGrid()
elapsed = (time.clock() - start)
logging.info("H5parm getValuesGrid() -- "+str(elapsed)+" s.")

print "Equal?", (Prot == Hrot3).all()

start = time.clock()
for val, axes in H.getIterValuesGrid():
    pass
elapsed = (time.clock() - start)
logging.info("H5parm getIterValuesGrid() -- "+str(elapsed)+" s.")

######################################################
# read whole file
logging.info("### Read and tabulate the whole file")
start = time.clock()
val, axes = H.getValuesGrid('')
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

######################################################
# read+write
logging.info("### Read all rotations for a dir/station and write them back")
Hw = solWriter(H5.getSoltab(solset,'amplitude000'))

start = time.clock()
Prot = P.getValues('CommonRotationAngle:CS001LBA')['CommonRotationAngle:CS001LBA']['values']
# parmdb write?
elapsed = (time.clock() - start)
logging.info("PARMDB -- "+str(elapsed)+" s.")

start = time.clock()
Hrot, axes, nrows = H.getValuesGrid(return_nrows = True)
Hw.setValueGrid(Hrot, nrows)
elapsed = (time.clock() - start)
logging.info("H5parm -- "+str(elapsed)+" s.")

del H
del P
