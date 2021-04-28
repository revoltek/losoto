#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os
import codecs
import numpy as np
from itertools import chain
from losoto.h5parm import h5parm, Soltab
from losoto import _version, _logging

_author = "Francesco de Gasperin (astro@voo.it)"

# get input
import argparse
parser = argparse.ArgumentParser(description='Combine h5parm files - '+_author)#, version=_version.__version__
parser.add_argument('h5parmFiles', nargs='+', help='List of h5parms')
parser.add_argument('--insolset', '-s', default='sol000', dest='insolset', help='Input solset name [default: sol000]')
parser.add_argument('--outsolset', '-u', default='sol000', dest='outsolset', help='Output solset name [default: sol000]')
parser.add_argument('--insoltab', '-t', default=None, dest='insoltab', help='Input soltab name (e.g. tabin000) or comma-separated list of soltab names - if not given use all')
parser.add_argument('--outh5parm', '-o', default='output.h5', dest='outh5parm', help='Output h5parm name [default: output.h5]')
parser.add_argument('--verbose', '-V', '-v', default=False, action='store_true', help='Go Vebose! (default=False)')
parser.add_argument('--squeeze', '-q', default=False, action='store_true', help='Remove all axes with a length of 1 (default=False)')
parser.add_argument('--history', '-H', default=False, action='store_true', help='Keep history of the export soltabs (default=False)')
parser.add_argument('--clobber', '-c', default=False, action='store_true', help='Replace exising outh5parm file instead of appending to it (default=False)')
args = parser.parse_args()

if len(args.h5parmFiles) < 1:
    parser.print_help()
    sys.exit(1)

# set logs
logger = _logging.Logger('info')
logging = _logging.logger
if args.verbose: logger.set_level("debug")

################################

# check input
if len(args.h5parmFiles) == 1 and (',' in args.h5parmFiles[0] or
                                   ('[' in args.h5parmFiles[0] and
                                    ']' in args.h5parmFiles[0])):
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
    if len(insoltabs) == 0:
        logging.critical('No soltabs found.')
        sys.exit()
else:
    insoltabs = args.insoltab.split(',')

# open input
h5s = []
for h5parmFile in args.h5parmFiles:
    h5 = h5parm(h5parmFile.replace("'",""), readonly=True)
    h5s.append(h5)

# open output
if os.path.exists(args.outh5parm) and args.clobber:
    os.remove(args.outh5parm)
h5Out = h5parm(args.outh5parm, readonly = False)

for insoltab in insoltabs:
    soltabs = []
    history = ''
    pointingNames = []; antennaNames = []
    pointingDirections = []; antennaPositions = []

    for h5 in h5s:
        solset = h5.getSolset(insolset)
        soltab = solset.getSoltab(insoltab)
        soltabs.append( soltab )
        history += soltab.getHistory()
        
        # collect pointings
        sous = solset.getSou()
        for k,v in list(sous.items()):
            if k not in pointingNames:
                pointingNames.append(k)
                pointingDirections.append(v)

        # collect anntennas
        ants = solset.getAnt()
        for k, v in list(ants.items()):
            if k not in antennaNames:
                antennaNames.append(k)
                antennaPositions.append(v)
                
    # create output axes
    logging.info("Sorting output axes...")
    axes = soltabs[0].getAxesNames()
    if args.squeeze:
        axes = [axis for axis in axes if soltabs[0].getAxisLen(axis) > 1 or axis == 'freq' ]
        removed_axes = list(set(soltabs[0].getAxesNames()) - set(axes))
        if len(removed_axes) == 0:
            args.squeeze = False
        else:
            axes_squeeze = tuple([soltabs[0].getAxesNames().index(removed_axis) for removed_axis in removed_axes ])
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
        if args.squeeze:
            allVals[np.ix_(*coords)] = np.squeeze(np.array(soltab.obj.val), axis = axes_squeeze)
            allWeights[np.ix_(*coords)] = np.squeeze(np.array(soltab.obj.weight), axis = axes_squeeze)
        else:
            allVals[np.ix_(*coords)] = soltab.obj.val
            allWeights[np.ix_(*coords)] = soltab.obj.weight
            

    # TODO: leave correct weights - this is a workaround for h5parm with weight not in float16
    allWeights[ allWeights != 0 ] = 1.

    # TODO: flag nans waiting for DPPP to do it
    allWeights[ np.isnan(allVals) ] = 0.

    logging.info('Writing output...')
    solsetOutName = args.outsolset
    soltabOutName = insoltab

    # create solset (and add all antennas and directions of other solsets)
    if solsetOutName in h5Out.getSolsetNames():
        solsetOut = h5Out.getSolset(solsetOutName)
    else:
        solsetOut = h5Out.makeSolset(solsetOutName)

    # create soltab
    soltabOut = solsetOut.makeSoltab(typ, soltabOutName, axesNames=axes, \
                     axesVals=[ allAxesVals[axis] for axis in axes ], \
                     vals=allVals, weights=allWeights)
    
    # add history table if requested
    if args.history:
        soltabOut.addHistory(history, date = False)

sourceTable = solsetOut.obj._f_get_child('source')
antennaTable = solsetOut.obj._f_get_child('antenna')
antennaTable.append(list(zip(*(antennaNames,antennaPositions))))
sourceTable.append(list(zip(*(pointingNames,pointingDirections))))

for h5 in h5s: h5.close()
logging.info(str(h5Out))
h5Out.close()
