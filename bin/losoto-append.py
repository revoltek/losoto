#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os
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
parser.add_argument('--intable', '-i', dest='intable', help='Input soltab name (e.g. sol000/tabin000)')
parser.add_argument('--outh5parm', '-o', dest='outh5parm', help='Output h5parm name')
parser.add_argument('--outtable', '-t', dest='outtable', help='Output soltab name (e.g. sol000/tabout000)')
parser.add_argument('--verbouse', '-v', dest='verbouse', action='store_true', help='Go Vebose! (default=False)')
parser.add_argument('--concataxis', '-c', dest='concataxis', help='Axis to concatenate (e.g. time)')
args = parser.parse_args()

if len(args.h5parms) < 2:
    argparse.print_help()
    sys.exit()

if args.verbose: _logging.setLevel("debug")

################################

class sfr(solFetcher):
    def __init__(self, table, concataxis)
        solFetcher.__init__(self, table = table, useCache = True)
        
        self.axesVals = []
        for axis in self.getAxesNames()
            self.axesVals.append( self.getAxisVals(axis) )

    def resample(self, axesValsNew):
        if axesValsNew != self.axesVals:
            from scipy.interpolate import griddata
            self.cacheVal = griddata(self.axesVals, self.cacheVal, axesValsNew, method='nearest')
            self.cacheWeight = griddata(self.axesVals, self.cacheWeight, axesValsNew, method='nearest')
            self.axesVals = axesValsNew

# read all tables
insolset, insoltab = args.intable.split('/')
sfs = []
pointingNames = []; antennaNames = []
pointingDirections = []; antennaPositions = []
for h5parmFile in args.h5parmFiles:
    logging.info("Reading "+h5parmFile)
    h5 = h5parm(h5parmFile, readonly = True)
    soltab = h5.getSoltab(solset=insolset, soltab=insoltab):
    sfs.append( sfr(soltab) )
    # collect pointings
    sous = h5.getSou(insolset)
    pointingNames.append(sous.keys())
    pointingDirections.append(sous.values())
    # collect anntennas
    ants = h5.getSou(insolset)
    antennaNames.append(ants.keys())
    antennaPositions.append(ants.values())

antennaNames = list(set(antennaNames))
antennaPositions = list(set(antennaPositions))
assert len(antennaNames) == len(antennaPositions) # different poistion for same antennas?!

# combine tables
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
        times.append(sf.getAxisVals('time'))
    if 'freq' in sf.getAxesNames():
        freqs.append(sf.getAxisVals('freq'))

axes = sfs[0].getAxesNames()
typ = sfs[0].getType()

if not args.concataxis in axes:
    logging.critical('Axis %s not found.' % args.concataxis)
    sys.exit(1)

# resampled time/freq axes values
if times != []:
    timeResamp = list(set(chain(*times)))
if freqs != []:
    freqResamp = list(set(chain(*freqs)))

# resampling of time/freq values
for sf in sfs:
    axesValsNew = sf.axesVals
    if 'time' in sf.getAxesNames():
        timeidx = axes.index('time')
        axesValsNew[timeidx] = timeResamp
    if 'freq' in sf.getAxesNames():
        freqidx = axes.index('freq')
        axesValsNew[freqidx] = freqResamp
    sf.resample(axesValsNew)

# sort tables on the concataxis first value
firstValsConcatAxis = [sf.getAxisValues(args.concataxis)[0] for sf in sfs]
sfs = [sf for (v,sf) in sorted(zip( firstValsConcatAxis, sfs))]

# re-order dirs (if concatenating in directions)
if args.concataxis == 'dir':
    pointingNames = [p for (v,p) in sorted(zip( firstValsConcatAxis, pointingNames))]
    pointingDirections = [p for (v,p) in sorted(zip( firstValsConcatAxis, pointingDirections))]

# get all vals/weights to concatenate
valsAll = []; weightsAll = []
for sf in sfs:
    valsAll.append( sf.getValues(retAxesVals=False) )
    weightsAll.append( sf.getValues(weights=True, retAxesVals=False) )

# creating concatenated table
axesVals = []; vals = []; weights = []
for a, axis in enumerate(axes):
    axisValsAll = []
    for i, sf in enumerate(sfs):
        val = sf.getAxisVals(axis)
        # allow concatenating directions all named 'pointing'
        if val == 'pointing':
            axisValsAll.append( 'pointing-%03i' % i )
            pointingNames[i] = 'pointing-%03i' % i
        else:
            axisValsAll.append( val )

    # check if all elements of the axis are equal
    if axisValsAll[1:] == axisValsAll[:-1]:
        axesVals.append(axisValsAll[0])

    elif axis == args.concataxis:
        # check intersection is empty
        if len( set.intersection(*map(set,axisValsAll)) ) > 0:
            logging.critical('Intersection of axis values in %s is not empty.' % axis)
        axesVals.append( [ [item for sublist in axisValsAll for item in sublist] ]) # flatten list of lists
        vals = np.concatenate(valsAll, axis=a)
        weights = np.concatenate(weightsAll, axis=a)

    else:
        # TODO: check missing axis values (e.g. missing antenna) and add flagged data
        logging.critical('Axis %s is not the same for all h5parm.' % axis)
        sys.exit(1)

# output
outsolset, outsoltab = args.outtable.split('/')
h5Out = h5parm(args.outh5parm, readonly = False)

# create solset (and add all antennas and directions of other solsets)
h5Out.makeSolset(outsolset)
sourceTable = solset._f_get_child('source')
antennaTable = solset._f_get_child('antenna')

# create soltab
h5Out.makeSoltab(outsolset, typ, outsoltab, axesNames=axes, \
                 axesVals=axesVals, vals=vals, weights=weights, parmdbType=sfs[0]._v_attrs['parmdb_type'])
sourceTable.append([(pointingNames,pointingDirections)])
antennaTable.append(zip(*(antennaNames,antennaPositions)))

del h5Out
