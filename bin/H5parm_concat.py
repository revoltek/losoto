#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, logging
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
args = parser.parse_args()

if len(args.h5parmFiles) < 2:
    argparse.print_help()
    sys.exit()

if args.verbose: _logging.setLevel("debug")

################################

class Soltabr(Soltab):
    def __init__(self, soltab, useCache = False, args = {}):
        Soltab.__init__(self, soltab, useCache = False, args = {})
        
    def resample(self, axisValsNew, axisName):
        """
        axisValsNew : new sampling values of the axis
        axisName : name of the axis to resample
        Get a list of axis values for an axisName and resample with neearest interpolation the values of the table
        """
        logging.info('Resampling...')
        axisIdx = self.getAxesNames().index(axisName)
        from scipy.interpolate import interp1d
        self.obj.val = interp1d(self.getAxisValues(axisName), self.obj.val, axis=axisIdx, kind='nearest', fill_value='extrapolate')(axisValsNew)
        self.obj.weight = interp1d(self.getAxisValues(axisName), self.obj.weight, axis=axisIdx, kind='nearest', fill_value='extrapolate')(axisValsNew)
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
insolset = args.insolset
if args.insoltab is None:
    # use all soltabs
    h5 = h5parm(args.h5parmFiles[0], readonly=True)
    solset = h5.getSolset(insolset)
    insoltabs = solset.getSoltabNames()
    h5.close()
else:
    insoltabs = [args.insoltab]

h5Out = h5parm(args.outh5parm, readonly = False)
for insoltab in insoltabs:
    soltabs = []
    pointingNames = []; antennaNames = []
    pointingDirections = []; antennaPositions = []
    for h5parmFile in args.h5parmFiles:
        logging.info("Reading "+h5parmFile)
        h5 = h5parm(h5parmFile, readonly=True)
        solset = h5.getSolset(insolset)
        soltab = Soltabr(solset.obj._f_get_child(insoltab)) # use inherited class 
        soltabs.append( soltab )
        # collect pointings
        sous = solset.getSou()
        pointingNames.append(sous.keys())
        pointingDirections.append(sous.values())
        # collect anntennas
        ants = solset.getAnt()
        [antennaNames.append(k) for k in ants.keys() if not k in antennaNames]
        [antennaPositions.append(v) for v in ants.values() if not any((v == x).all() for x in antennaPositions)]
    
    # combine tables
    logging.info('Ordering data...')
    times = []
    freqs = []
    for soltab in soltabs:
        if soltab.getAxesNames() != soltabs[0].getAxesNames():
            logging.critical('Input soltabs have different axes.')
            sys.exit(1)
        if soltab.getType() != soltabs[0].getType():
            logging.critical('Input soltabs have different types.')
            sys.exit(1)
        if 'time' in soltab.getAxesNames():
            times.append(soltab.getAxisValues('time'))
        if 'freq' in soltab.getAxesNames():
            freqs.append(soltab.getAxisValues('freq'))
    
    axes = soltabs[0].getAxesNames()
    typ = soltabs[0].getType()
    
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
        print '- Will be:', len(timeResamp)
    if freqs != []:
        freqResamp = np.array(sorted(list(set(chain(*freqs)))))
        print 'len freqs:',
        for f in freqs:
            print len(f),
        print '- Will be:', len(freqResamp)
    
    # resampling of time/freq values
    for soltab in soltabs:
        if 'time' in axes and len(timeResamp) != len(soltab.getAxisValues('time')) and args.concataxis != 'time':
            soltab.resample(timeResamp, 'time')
        # first resample in time, then in freq
        if 'freq' in axes and len(freqResamp) != len(soltab.getAxisValues('freq')) and args.concataxis != 'freq':
            soltab.resample(freqResamp, 'freq')
    
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
            logging.critical('Axis %s is not the same for all h5parm.' % axis)
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
