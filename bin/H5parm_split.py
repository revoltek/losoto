#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This tool is used to split an H5parm into multiple smaller ones.

import sys, os
import numpy as np
from losoto.h5parm import h5parm
from losoto import _version, _logging

_author = "David Rafferty (drafferty@hs.uni-hamburg.de)"

# get input
import argparse
parser = argparse.ArgumentParser(description='Split h5parm files - '+_author)#, version=_version.__version__
parser.add_argument('h5parmFile', help='Input h5parm')
parser.add_argument('axis', help='Axis name to split')
parser.add_argument('length', help='Target length of axis in each output h5parm (in the same units as the axis)')
parser.add_argument('--insolset', '-s', default='sol000', dest='insolset', help='Input solset name [default: sol000]')
parser.add_argument('--outsolset', '-o', default=None, dest='outsolset', help='Output solset name [default: same as insolset]')
parser.add_argument('--root', '-r', default='output', dest='root', help='Root of output h5parm names; e.g., root_0.h5, root_1.h5, etc. [default: output]')
parser.add_argument('--verbose', '-V', '-v', default=False, action='store_true', help='Go Vebose! (default=False)')
parser.add_argument('--clobber', '-c', default=False, action='store_true', help='Replace exising outh5parm files (default=False)')
args = parser.parse_args()

if len(args.h5parmFile) < 1:
    parser.print_help()
    sys.exit(1)

# log
logger = _logging.Logger('info')
logging = _logging.logger
if args.verbose: logger.set_level("debug")

################################

# read input table
h5 = h5parm(args.h5parmFile, readonly=True)
solset = h5.getSolset(args.insolset)
insoltabs = solset.getSoltabNames()
soltabs = []
for insoltab in insoltabs:
    soltab = solset.getSoltab(insoltab)
    soltabs.append( soltab )
pointingNames = []; antennaNames = []
pointingDirections = []; antennaPositions = []
ants = solset.getAnt()
sous = solset.getSou()
for k,v in list(sous.items()):
    if k not in pointingNames:
        pointingNames.append(k)
        pointingDirections.append(v)
for k, v in list(ants.items()):
    if k not in antennaNames:
        antennaNames.append(k)
        antennaPositions.append(v)

# check that all soltabs have the requested axis and decide how many output files to make
for s in soltabs:
    if args.axis not in s.getAxesNames():
        logging.critical("Axis '{0}' not present in soltab '{1}'.".format(args.axis, s.name))
        sys.exit(1)

    vals = s.getAxisValues(args.axis)

    # delineate the chunks
    chunk_ind = np.where(((vals-vals[0]) % float(args.length))[1:] < np.diff(vals))[0]
    nchunks = len(chunk_ind) + 1

# fill the output files
for i in range(nchunks):
    outh5file = '{0}_{1}.h5'.format(args.root, i)
    if os.path.exists(outh5file) and args.clobber:
        os.remove(outh5file)
    logging.info("\nCreating h5parm {}".format(outh5file))
    outh5 = h5parm(outh5file, readonly=False)
    solsetOutName = args.outsolset
    if solsetOutName is None:
        solsetOutName = args.insolset
    solsetOut = outh5.makeSolset(solsetOutName)

    for s in soltabs:
        vals = s.getAxisValues(args.axis)
        if i == 0:
            startval = vals[0]
        else:
            startval = vals[chunk_ind[i-1]+1]
        if i == nchunks-1:
            endval = vals[-1]
        else:
            endval = vals[chunk_ind[i]+1]
        s.setSelection(**{args.axis:{'min':startval,'max':endval,'step':1}})
        solsetOut.makeSoltab(s.getType(), s.name, axesNames=s.getAxesNames(),
                         axesVals=[s.getAxisValues(a) for a in s.getAxesNames()],
                         vals=s.getValues(retAxesVals=False),
                         weights=s.getValues(weight=True, retAxesVals=False))
        s.clearSelection()

    sourceTable = solsetOut.obj._f_get_child('source')
    antennaTable = solsetOut.obj._f_get_child('antenna')
    antennaTable.append(list(zip(*(antennaNames,antennaPositions))))
    sourceTable.append(list(zip(*(pointingNames,pointingDirections))))
    outh5.close()
