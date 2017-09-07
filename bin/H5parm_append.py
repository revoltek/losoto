#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Authors:
# Francesco de Gasperin
# David Raffery
# Maaijke Mevius
_author = "Francesco de Gasperin (fdg@strw.leidenuniv.nl)"

import sys, os, logging
import numpy as np
from itertools import chain
from losoto.h5parm import solFetcher
from losoto.h5parm import h5parm
from losoto import _version
from losoto import _logging

# get input
import argparse
parser = argparse.ArgumentParser(description='Combine tables in h5parm files')
parser.add_argument('h5parmFiles', nargs='+', help='List of h5parms')
parser.add_argument('--insolsset', '-s', default='sol000', dest='insolset', help='Input solset name [default: sol000]')
parser.add_argument('--insoltab', '-t', dest='insoltab', help='Input soltab name (e.g. tabin000) - if not give use all')
parser.add_argument('--outh5parm', '-o', dest='outh5parm', help='Output h5parm name')
parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', help='Go Vebose! (default=False)')
parser.add_argument('--concataxis', '-c', dest='concataxis', help='Axis to concatenate (e.g. time)')
args = parser.parse_args()

if len(args.h5parmFiles) < 2:
    argparse.print_help()
    sys.exit()

if args.verbose: _logging.setLevel("debug")

################################

class sfr(solFetcher):
    def __init__(self, table):
        solFetcher.__init__(self, table = table, useCache = True)
        
    def resample(self, axisValsNew, axisName):
        """
        axisValsNew : new sampling values of the axis
        axisName : name of the axis to resample
        Get a list of axis values for an axisName and resample with neearest interpolation the values of the table
        """
        logging.info('Resampling...')
        axisIdx = self.getAxesNames().index(axisName)
        from scipy.interpolate import interp1d
        self.cacheVal = interp1d(self.getAxisValues(axisName), self.cacheVal, axis=axisIdx, kind='nearest', fill_value='extrapolate')(axisValsNew)
        self.cacheWeight = interp1d(self.getAxisValues(axisName), self.cacheWeight, axis=axisIdx, kind='nearest', fill_value='extrapolate')(axisValsNew)
        # update axis vals
        self.axes[axisName] = axisValsNew
        self.setSelection() # reset selection to empty since values are different

def equalArr(arr):
    """
    Check if arrays are equal
    """
    for a in arr[1:]:
        if not np.array_equal(a, arr[0]):
            return False
    return True

# read all tables
insolset = args.insolsset
insoltab = args.insoltab
sfs = []
pointingNames = []; antennaNames = []
pointingDirections = []; antennaPositions = []
for h5parmFile in args.h5parmFiles:
    logging.info("Reading "+h5parmFile)
    h5 = h5parm(h5parmFile, readonly = True)
    soltab = h5.getSoltab(solset=insolset, soltab=insoltab)
    sfs.append( sfr(soltab) )
    # collect pointings
    sous = h5.getSou(insolset)
    pointingNames.append(sous.keys())
    pointingDirections.append(sous.values()[0])
    # collect anntennas
    ants = h5.getAnt(insolset)
    [antennaNames.append(k) for k in ants.keys() if not k in antennaNames]
    [antennaPositions.append(v) for v in ants.values() if not any((v == x).all() for x in antennaPositions)]

# combine tables
logging.info('Ordering data...')
times = []
freqs = []
for sf in sfs:
    if sf.getAxesNames() != sfs[0].getAxesNames():
        logging.critical('Input soltabs have different axes.')
        sys.exit(1)
    if sf.getType() != sfs[0].getType():
        logging.critical('Input soltabs have different types.')
        sys.exit(1)
    if 'time' in sf.getAxesNames():
        times.append(sf.getAxisValues('time'))
    if 'freq' in sf.getAxesNames():
        freqs.append(sf.getAxisValues('freq'))

axes = sfs[0].getAxesNames()
typ = sfs[0].getType()

if not args.concataxis in axes:
    logging.critical('Axis %s not found.' % args.concataxis)
    sys.exit(1)

# resampled time/freq axes values
# every single time/freq valu for all tables is in these arrays (ordered)
if times != []:
    timeResamp = np.array(sorted(list(set(chain(*times)))))
    print 'len times:',
    for t in times:
        print len(t),
    print 'Resamp to:', len(timeResamp)
if freqs != []:
    freqResamp = np.array(sorted(list(set(chain(*freqs)))))
    print 'len freqs:',
    for f in freqs:
        print len(f),
    print 'Resamp to:', len(freqResamp)

# resampling of time/freq values
for sf in sfs:
    if 'time' in sf.getAxesNames() and len(timeResamp) != len(sf.getAxisValues('time')):
        sf.resample(timeResamp, 'time')
    # first resample in time, then in freq
    if 'freq' in sf.getAxesNames() and len(freqResamp) != len(sf.getAxisValues('freq')):
        sf.resample(freqResamp, 'freq')

# sort tables based on the first value of the concatAxis
logging.info('Sorting tables...')
firstValsConcatAxis = [sf.getAxisValues(args.concataxis)[0] for sf in sfs]
idxToSort = [i[0] for i in sorted(enumerate(firstValsConcatAxis), key=lambda x:x[1])]
print 'idxToSort', idxToSort
sfs = [sfs[i] for i in idxToSort]

# re-order dirs (if concatenating in directions)
if args.concataxis == 'dir':
    pointingNames = [pointingNames[i] for i in idxToSort]
    pointingDirections = [pointingDirections[i] for i in idxToSort]

# get all vals/weights to concatenate
logging.info('Allocate final data...')
valsAll = []; weightsAll = []
for sf in sfs:
    valsAll.append( sf.cacheVal )
    weightsAll.append( sf.cacheWeight )

# creating concatenated table
logging.info('Concatenate final dataset...')
axesVals = []; vals = []; weights = []
for a, axis in enumerate(axes):
    axisValsAll = []
    for i, sf in enumerate(sfs):
        val = sf.getAxisValues(axis)
        # allow concatenating directions all named 'pointing'
        if axis == 'dir' and len(val) == 1 and val[0] == 'pointing':
            axisValsAll.append( np.array(['pointing-%03i' % i]) )
            pointingNames[i] = 'pointing-%03i' % i
        else:
            axisValsAll.append( val )

    # check if all elements of the axis are equal
    if equalArr(axisValsAll):
        axesVals.append(axisValsAll[0])

    elif axis == args.concataxis:
        # check intersection is empty
        if len( set.intersection(*map(set,axisValsAll)) ) > 0:
            logging.critical('Intersection of axis values in %s is not empty.' % axis)
        axesVals.append( np.array( [item for sublist in axisValsAll for item in sublist] ) ) # flatten list of lists
        vals = np.concatenate(valsAll, axis=a)
        weights = np.concatenate(weightsAll, axis=a)

    else:
        # TODO: check missing axis values (e.g. missing antenna) and add flagged data
        logging.critical('Axis %s is not the same for all h5parm.' % axis)
        sys.exit(1)

logging.info('Writing output...')
# output
outsolset = args.insolset
outsoltab = args.insoltab
h5Out = h5parm(args.outh5parm, readonly = False)

# create solset (and add all antennas and directions of other solsets)
solset = h5Out.makeSolset(outsolset)
sourceTable = solset._f_get_child('source')
antennaTable = solset._f_get_child('antenna')

# create soltab
h5Out.makeSoltab(outsolset, typ, outsoltab, axesNames=axes, \
                 axesVals=axesVals, vals=vals, weights=weights, parmdbType=soltab._v_attrs['parmdb_type'])

antennaTable.append(zip(*(antennaNames,antennaPositions)))
sourceTable.append(zip(*(pointingNames,pointingDirections)))

del h5Out
