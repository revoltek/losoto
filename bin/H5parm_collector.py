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
args = parser.parse_args()

if len(args.h5parmFiles) < 2:
    parser.print_help()
    sys.exit()

if args.verbose: _logging.setLevel("debug")

################################

# check input
if len(args.h5parmFiles) == 1 and ',' in args.h5parmFiles[0]:
    # Treat input as a string with comma-separated values
    args.h5parmFiles = args.h5parmFiles[0].strip('[]').split(',')
    args.h5parmFiles  = [f.strip() for f in args.h5parmFiles]

# read all tables
insolset = args.insolset
if args.insoltab is None:
    # use all soltabs, find out names from first h5parm
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
        soltab = solset.getSoltab(insoltab)
        soltabs.append( soltab )
        # collect pointings
        sous = solset.getSou()
        pointingNames.extend(sous.keys())
        pointingDirections.extend(sous.values())
        # collect anntennas
        ants = solset.getAnt()
        [antennaNames.append(k) for k in ants.keys() if not k in antennaNames]
        [antennaPositions.append(v) for v in ants.values() if not any((v == x).all() for x in antennaPositions)]

    # create output axes
    logging.info("Sorting output axes...")
    axes = soltabs[0].getAxesNames()
    typ = soltabs[0].getType()
    allAxesVals = {axis:[] for axis in axes}
    allShape = []
    for axis in axes:
        print('len %s:' % axis, end='')
        for soltab in soltabs:
            allAxesVals[axis].append( soltab.getAxisValues(axis) )
            print(' %i' % soltab.getAxisLen(axis), end='')
        allAxesVals[axis] = np.array(sorted(list(set(chain(*allAxesVals[axis])))))
        allShape.append(len(allAxesVals[axis]))
        print(' - Will be: %i' % len(allAxesVals[axis]))

    # make final arrays
    logging.info("Allocating space...")
    logging.debug("Shape:"+str(allShape))
    allVals = np.empty( shape=allShape )
    allVals[:] = np.nan
    allWeights = np.zeros( shape=allShape )#, dtype=np.float16 )

    # fill arrays
    logging.info("Filling new table...")
    for soltab in soltabs:
        coords = []
        for axis in axes:
            coords.append( np.searchsorted( allAxesVals[axis], soltab.getAxisValues(axis) ) )
        allVals[np.ix_(*coords)] = soltab.obj.val
        allWeights[np.ix_(*coords)] = soltab.obj.weight

    # TODO: leave correct weights - this is a workaround for h5parm with weight not in float16
    allWeights[ allWeights != 0 ] = 1.

    # TODO: interpolate on selected axes

    logging.info('Writing output...')
    solsetOutName = args.insolset
    soltabOutName = args.insoltab

    # create solset (and add all antennas and directions of other solsets)
    if solsetOutName in h5Out.getSolsetNames():
        solsetOut = h5Out.getSolset(solsetOutName)
    else:
        solsetOut = h5Out.makeSolset(solsetOutName)

    # create soltab
    solsetOut.makeSoltab(typ, soltabOutName, axesNames=axes, \
                     axesVals=[ allAxesVals[axis] for axis in axes ], \
                     vals=allVals, weights=allWeights)

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
