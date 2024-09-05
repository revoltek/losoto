#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
from scipy.interpolate import interp1d
from losoto.h5parm import h5parm, Soltab
from losoto import _version, _logging

_author = "Francesco de Gasperin (astro@voo.it)"

# get input
import argparse
parser = argparse.ArgumentParser(description='Interpolate h5parm files on the fastest varying - '+_author)#, version=_version.__version__
parser.add_argument('h5parmFiles', nargs='+', help='List of h5parms')
parser.add_argument('--insolset', '-s', default='sol000', dest='insolset', help='Input solset name [default: sol000]')
parser.add_argument('--outsolset', '-u', default='sol000', dest='outsolset', help='Output solset name [default: sol000]')
parser.add_argument('--insoltab', '-t', default=None, dest='insoltab', help='Input soltab name (e.g. tabin000) or comma-separated list of soltab names - if not given use all')
parser.add_argument('--interpaxes', '-a', default='time,freq', dest='interpaxes', help='Comma-separated list of axes to interpolate, all other axes need to be equal - default: time,freq')
parser.add_argument('--outh5parm', '-o', default='output.h5', dest='outh5parm', help='Output h5parm name [default: output.h5]')
parser.add_argument('--method', '-m', default='nearest', dest='method', help='see scipy interp1d man, possibilities: linear, nearest, zero, slinear, quadratic, cubic, previous, next. Default: nearest.')
parser.add_argument('--verbose', '-V', '-v', default=False, action='store_true', help='Go Vebose! (default=False)')
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

# open input
all_soltabs = {} # {'amplitude':[soltab1, soltab2], 'phase':[soltab3]}
all_axes_names = {} # {'amplitude':['dir','time','ant','freq'], 'phase':['dir','time','freq']} # dir is the first by construction
all_axes_vals = {} # {'time':[[1,2,34,5],[1,3,5]],'freq':[[1,2,3,4,5],[1,3,5]]}
pointingNames = []; antennaNames = []; pointingDirections = []; antennaPositions = []
for h5parmFile in args.h5parmFiles:
    logging.info('Working on %s...' % h5parmFile)
    h5 = h5parm(h5parmFile.replace("'",""), readonly=True)
    solset = h5.getSolset(args.insolset)
    insoltabs = solset.getSoltabNames()
    # restrict to requested insoltabs
    if args.insoltab is not None:
        insoltabs = [insoltab for insoltab in insoltabs if insoltab in args.insoltab]
    if len(insoltabs) == 0:
        logging.critical('No soltabs found.')
        sys.exit()
    for insoltab in insoltabs:
        soltab = solset.getSoltab(insoltab)

        # classify soltabs
        typ = soltab.getType()
        if typ in all_soltabs.keys():
            all_soltabs[typ].append( soltab )
        else:
            all_soltabs[typ] = [soltab]

        # classify axes
        axes = soltab.getAxesNames()
        for axis in axes:
            if axis in all_axes_vals.keys():
                all_axes_vals[axis].append( soltab.getAxisValues(axis) )
            else:
                all_axes_vals[axis] = [soltab.getAxisValues(axis)]

        # add axes names for this table typ
        axes = ['dir']+[axis for axis in axes if axis != 'dir']
        if typ in all_axes_names.keys():
            if all_axes_names[typ] != axes:
                logging.error('Table %s has different axes from previously loadded tables of the same type.' % insoltab)
                sys.exit()
        else:
            all_axes_names[typ] = axes


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
   
# Find fastest varying (longest) axis for each axis type
fast_axes_vals = {}
for axis, all_vals in all_axes_vals.items():
    if axis == 'dir':
        fast_axes_vals['dir'] = list(set([val[0] for val in all_vals]))
    else:
        logging.info('Looking for fastest axis in %s...' % axis)
        lens = [len(vals) for vals in all_vals]
        fast_axes_vals[axis] = all_vals[lens.index(max(lens))]
        print(lens,'->',len(fast_axes_vals[axis]))
        if axis not in ['time','freq'] and any([len(fast_axes_vals[axis]) != x for x in lens]):
            logging.error('Only time and freq can be resampled, not %s.' % axis)
            sys.exit()
 
# open output
if os.path.exists(args.outh5parm) and args.clobber:
    os.remove(args.outh5parm)
h5Out = h5parm(args.outh5parm, readonly = False)
# create solset (and add all antennas and directions of other solsets)
solsetOut = h5Out.makeSolset(args.outsolset)
sourceTable = solsetOut.obj._f_get_child('source')
antennaTable = solsetOut.obj._f_get_child('antenna')
antennaTable.append(list(zip(*(antennaNames,antennaPositions))))
sourceTable.append(list(zip(*(pointingNames,pointingDirections))))

for typ, soltabs in all_soltabs.items():
    logging.info('Creating new table of type: %s...' % typ)
    shape = [len(fast_axes_vals[axis]) for axis in all_axes_names[typ]]
    print('Shape: ', shape, ' - ', all_axes_names[typ])
    if typ == 'amplitude':
        new_vals = np.ones(shape=shape)
    else:
        new_vals = np.zeros(shape=shape)
    new_weights = np.ones(shape=shape)

    for soltab in soltabs:
        vals = soltab.getValues(retAxesVals=False)
        weights = soltab.getValues(retAxesVals=False, weight=True)
        #vals[np.isnan(vals)] = 0 # avoid nans (these are flagged data)
        for axis_idx, axis in enumerate(soltab.getAxesNames()):
            if axis not in ['time','freq']: continue
            # get axies
            axis_vals = soltab.getAxisValues(axis)
            new_axis_vals = fast_axes_vals[axis]

            if len(axis_vals) < len(new_axis_vals):
                if len(axis_vals) == 1:
                    pass # numpy should do the proper casting
                else:
                    # interpolate on faster varying axis
                    f = interp1d(axis_vals, vals, kind=args.method, axis=axis_idx, copy=True, fill_value='extrapolate', assume_sorted=True)
                    vals = f(new_axis_vals)
                    fw = interp1d(axis_vals, weights, kind=args.method, axis=axis_idx, copy=True, fill_value='extrapolate', assume_sorted=True)
                    weights = fw(new_axis_vals)

        # remove direciton axis
        idx_axis_dir = soltab.getAxesNames().index('dir') # position of the direction axis
        vals = np.squeeze(vals, axis=idx_axis_dir)
        weights = np.squeeze(weights, axis=idx_axis_dir)
        # find direction
        idx_dir = fast_axes_vals['dir'].index(soltab.getAxisValues('dir'))
        # add to global value
        if typ == 'amplitude':
            new_vals[idx_dir] *= vals
        else:
            new_vals[idx_dir] += vals
        new_weights[idx_dir] *= weights

    # write new soltab (default name)
    soltabOut = solsetOut.makeSoltab(typ, axesNames=all_axes_names[typ],
                axesVals=[ fast_axes_vals[axis_name] for axis_name in all_axes_names[typ] ],
                vals=new_vals, weights=new_weights)

logging.info(str(h5Out))
h5Out.close()
