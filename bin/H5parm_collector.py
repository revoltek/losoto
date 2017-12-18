#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys, os, logging
import codecs
import numpy as np
from itertools import chain
from losoto.h5parm import h5parm, Soltab
from losoto import _version, _logging

_author = "Francesco de Gasperin (astro@voo.it)"

# get input
import argparse
parser = argparse.ArgumentParser(description='Combine h5parm files - '+_author, version=_version.__version__)
parser.add_argument('h5parmFiles', nargs='+', help='List of h5parms')
parser.add_argument('--insolset', '-s', default='sol000', dest='insolset', help='Input solset name [default: sol000]')
parser.add_argument('--insoltab', '-t', default=None, dest='insoltab', help='Input soltab name (e.g. tabin000) - if not given use all')
parser.add_argument('--outh5parm', '-o', default='output.h5', dest='outh5parm', help='Output h5parm name [default: output.h5]')
parser.add_argument('--verbose', '-V', dest='verbose', action='store_true', help='Go Vebose! (default=False)')
parser.add_argument('--concataxis', '-c', dest='concataxis', help='Axis to concatenate (e.g. time)')
parser.add_argument('--resampaxes', '-r', nargs='*', dest='resampaxes', default=['time','freq'], help='Axes to resanple (default: time,freq)') # nargs=* means 0 or more
parser.add_argument('--fillaxes', '-f', nargs='*', dest='fillaxes', default=['ant'], help='Axes to fill with flagged data for missing values (default: ant)')
args = parser.parse_args()

if len(args.h5parmFiles) < 1:
    parser.print_help()
    sys.exit()

if args.verbose: _logging.setLevel("debug")

################################

class Soltabr(Soltab):
    def __init__(self, soltab, useCache = False, args = {}):
        Soltab.__init__(self, soltab, useCache = False, args = {})

    def resample(self, axisValsNew, axisName, flag=False):
        """
        Get a list of axis values for an axisName and resample with neearest interpolation the values of the table
        axisValsNew : new sampling values of the axis
        axisName : name of the axis to resample
        flag: if True the expanded values are flagged
        """
        logging.info('Resampling: %s...' % axisName)
        axisIdx = self.getAxesNames().index(axisName)
        from scipy.interpolate import interp1d

        # get axis values and update them
        axisVals = self.getAxisValues(axisName)
        self.axes[axisName] = axisValsNew

        # transform str axis into numbers, this is needed e.g. for antennas
        if isinstance(axisValsNew[0], str):
            axisValsNew = [int(codecs.encode(s, 'hex'), 16) for s in axisValsNew]
            axisVals = [int(codecs.encode(s, 'hex'), 16) for s in axisVals]

        self.obj.val = interp1d(axisVals, self.obj.val, axis=axisIdx, kind='nearest', fill_value='extrapolate')(axisValsNew)
        self.obj.weight = interp1d(axisVals, self.obj.weight, axis=axisIdx, kind='nearest', fill_value='extrapolate')(axisValsNew)

        # flag all resampled data
        if flag:
            # move relevant axis at the end
            vals = np.swapaxes(self.obj.val, axisIdx, 0)
            weights = np.swapaxes(self.obj.weight, axisIdx, 0)
            for i, val in enumerate(axisValsNew):
                if not val in axisVals:
                    weights[i] = 0.
                    vals[i] = np.nan
            self.obj.val = np.swapaxes(vals, axisIdx, 0)
            self.obj.weight = np.swapaxes(weights, axisIdx, 0)

        self.setSelection() # reset selection to empty since values are different

def equalArr(arr):
    """
    Check if arrays are equal
    """
    for a in arr[1:]:
        if not np.array_equal(a, arr[0]):
            return False
    return True

# check input
if len(args.h5parmFiles) == 1 and ',' in args.h5parmFiles[0]:
    # Treat input as a string with comma-separated values
    args.h5parmFiles = args.h5parmFiles[0].strip('[]').split(',')
    args.h5parmFiles  = [f.strip() for f in args.h5parmFiles]

if args.resampaxes is None: args.resampaxes = []
if args.fillaxes is None: args.fillaxes = []

# read all tables
insolset = args.insolset
if args.insoltab is None:
    # use all soltabs
    h5 = h5parm(args.h5parmFiles[0], readonly=True)
    solset = h5.getSolset(insolset)
    insoltabs = solset.getSoltabNames()
    h5.close()
else:
    insoltabs = [args.insoltab]

# open input
h5s = []
for h5parmFile in args.h5parmFiles:
    h5 = h5parm(h5parmFile, readonly=True)
    h5s.append(h5)
# open output
h5Out = h5parm(args.outh5parm, readonly = False)

for insoltab in insoltabs:
    soltabs = []
    pointingNames = []; antennaNames = []
    pointingDirections = []; antennaPositions = []

    for h5 in h5s:
        solset = h5.getSolset(insolset)
        soltab = Soltabr(solset.obj._f_get_child(insoltab)) # use inherited class
        soltabs.append( soltab )
        # collect pointings
        sous = solset.getSou()
        pointingNames.extend(sous.keys())
        pointingDirections.extend(sous.values())
        # collect anntennas
        ants = solset.getAnt()
        [antennaNames.append(k) for k in ants.keys() if not k in antennaNames]
        [antennaPositions.append(v) for v in ants.values() if not any((v == x).all() for x in antennaPositions)]

    # combine tables
    logging.info('Ordering data...')
    resamp = {axis:[] for axis in args.fillaxes+args.resampaxes} # dict to store+compute resamp/fill axis values
    for soltab in soltabs:
        if soltab.getAxesNames() != soltabs[0].getAxesNames():
            logging.critical('Input soltabs have different axes.')
            sys.exit(1)
        if soltab.getType() != soltabs[0].getType():
            logging.critical('Input soltabs have different types.')
            sys.exit(1)
        # concatenate axis values inside the resamp dict
        for axis in soltab.getAxesNames():
            if axis in resamp.keys():
                resamp[axis].append(soltab.getAxisValues(axis))

    axes = soltabs[0].getAxesNames()
    typ = soltabs[0].getType()

    if not args.concataxis in axes:
        logging.critical('Axis %s not found.' % args.concataxis)
        sys.exit(1)

    # all arrays in each entry of the resamp dict are combined in a single one
    # this array has all values, not repeated and sorted
    for axis, allAxisVals in iter(resamp.items()):
        resamp[axis] = np.array(sorted(list(set(chain(*allAxisVals)))))
        print('len %s:' % axis, end='')
        for axisVals in allAxisVals:
            print(' %i' % len(axisVals), end='')
        print('Will be: %i' % len(resamp[axis]))

    # resampling/filling of values
    for soltab in soltabs:
        for axis in soltab.getAxesNames():
            if axis in resamp.keys():
                if axis in args.resampaxes and len(resamp[axis]) != len(soltab.getAxisValues(axis)) and args.concataxis != axis:
                    soltab.resample(resamp[axis], axis, flag=False)
                elif axis in args.fillaxes and len(resamp[axis]) != len(soltab.getAxisValues(axis)) and args.concataxis != axis:
                    soltab.resample(resamp[axis], axis, flag=True)

    # sort tables based on the first value of the concatAxis
    logging.info('Sorting tables...')
    firstValsConcatAxis = [soltab.getAxisValues(args.concataxis)[0] for soltab in soltabs]
    idxToSort = [i[0] for i in sorted(enumerate(firstValsConcatAxis), key=lambda x:x[1])]
    #print 'idxToSort', idxToSort
    soltabs = [soltabs[i] for i in idxToSort]

    # re-order dirs (if concatenating in directions)
    if args.concataxis == 'dir':
        pointingNames = [pointingNames[i] for i in idxToSort]
        pointingDirections = [pointingDirections[i] for i in idxToSort]

    # get all vals/weights to concatenate
    logging.info('Allocate final data...')
    valsAll = []; weightsAll = []
    for soltab in soltabs:
        valsAll.append( soltab.obj.val )
        weightsAll.append( soltab.obj.weight )

    # creating concatenated table
    logging.info('Concatenate final dataset...')
    axesVals = []; vals = []; weights = []
    for a, axis in enumerate(axes):
        axisValsAll = []
        for i, soltab in enumerate(soltabs):
            val = soltab.getAxisValues(axis)
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
            logging.critical('Axis %s is not the same for all h5parms.' % axis)
            sys.exit(1)

    logging.info('Writing output...')
    solsetOutName = args.insolset
    soltabOutName = args.insoltab

    # create solset (and add all antennas and directions of other solsets)
    if solsetOutName in h5Out.getSolsetNames():
        solsetOut = h5Out.getSolset(solsetOutName)
    else:
        solsetOut = h5Out.makeSolset(solsetOutName)

    # create soltab
    weights[ np.isnan(vals) ] = 0. # to be sure, can be removed when DPPP does it properly
    logging.debug( "Set weight to 0 for %i values." % np.sum(np.isnan(vals)) )
    solsetOut.makeSoltab(typ, soltabOutName, axesNames=axes, \
                     axesVals=axesVals, vals=vals, weights=weights)

sourceTable = solsetOut.obj._f_get_child('source')
antennaTable = solsetOut.obj._f_get_child('antenna')
try:
    antennaTable.append(zip(*(antennaNames,antennaPositions)))
except:
    logging.warning('Couldnt fill antenna table.')
try:
    sourceTable.append(zip(*(pointingNames,pointingDirections)))
except:
    logging.warning('Couldnt fill source table.')

for h5 in h5s: h5.close()
logging.info(str(h5Out))
h5Out.close()
