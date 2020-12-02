#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation interpolates the directions of a soltab, adding more directions


import numpy as np
import scipy
from scipy.spatial.distance import cdist
from scipy.interpolate import Rbf
from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading INTERPOLATEDIRECTIONS module.')

# this funct is called by losoto to set parameters and call the real run()
def _run_parser(soltab, parser, step):
    interp_dirs = parser.getarrayfloat2d( step, 'interp_dirs') # no default
    soltabOut = parser.getstr( step, 'soltabOut', "")
    prefix = parser.getstr( step, 'prefix', 'interp_')
    ncpu = parser.getint( '_global', 'ncpu', 0 )

    parser.checkSpelling( step, soltab, ['interp_dirs', 'soltabOut', 'prefix', 'ncpu'])
    return run(soltab, interp_dirs, soltabOut, prefix, ncpu)


def _haversine(s1, s2):
    """
    Calculate the great circle distance between two points
    (specified in rad)
    """
    return 2*np.arcsin(np.sqrt(np.sin((s2[1]-s1[1])/2.0)**2 + np.cos(s1[1]) * np.cos(s2[1]) * np.sin((s2[0]-s1[0])/2.0)**2))

def _cartesian_from_radec(x):
    return np.array([np.cos(x[:,0])*np.sin(x[:,1]), np.sin(x[:,0])*np.sin(x[:,1]), np.cos(x[:,1])])

def interpolate_directions2d(cal_vals, cal_weights, cal_dirs, interp_dirs, interp_kind, smooth=0.0):
    """
    Function to parallelize interpolation. Interpolates and evaluates
    the interpolation at the new directions for one ant/pol/time etc.

    Parameters
    ----------
    cal_vals: array, values of selection, 1st axis: dir, other axes flattened
    weights: weights of values,
    cal_dirs: (2,n) array of floats, ra/dec values of calibrator directions
    interp_dirs: (2,n) array of floats, ra/dec values of interpolate directions
    interp_kind: string, 'wrap' or 'lin' dependent on wheter to interpolate phases
                 (independently in real and imag) or TEC, amp..

    Returns
    -------
    new_vals: array of floats, interpolated values of selection at new directions
    new_weights: array of floats, weights
    """
    # approximation: Use euclidean norm in 3d instead of haversine. Every implementation of
    # haversine I found is much slower than scipy.spatial.distance.cdist for eculidean
    # to 3d for euclidean approximation. 2d-interpolating in ra-dec has the issue of wraps.

    new_weights = np.ones((len(interp_dirs),cal_vals.shape[-1]))
    cal_weights[np.isnan(cal_vals)] = 0.0 # make sure all NaNs are flagged
    cal_vals[np.isnan(cal_vals)] = 0.0 # set flagged values to 0
    new_weights[:,np.sum(cal_weights, axis=0) < len(cal_weights)] # Set new weights to 0 where at least one old direction is flagged.

    # cal_dirs = _cartesian_from_radec(cal_dirs).T
    # interp_dirs = _cartesian_from_radec(interp_dirs).T
    if interp_kind == 'wrap': # for phases, interpolate real and imag.
        _complex = np.exp(1.j * cal_vals)
        f_interp_re = Rbf(*cal_dirs.T, _complex.real, norm=_haversine, mode='N-D', smooth=smooth) # multiquadratic, linear
        f_interp_im = Rbf(*cal_dirs.T, _complex.imag, norm=_haversine, mode='N-D', smooth=smooth)
        interp_re = f_interp_re(*interp_dirs.T)
        interp_im = f_interp_im(*interp_dirs.T)
        new_vals = np.angle(interp_re + 1.j * interp_im)
    elif interp_kind == 'lin':
        f_interp = Rbf(*cal_dirs.T, cal_vals, norm=_haversine, mode='N-D', smooth=smooth)
        new_vals = f_interp(*interp_dirs.T)
    return new_vals, new_weights

def interpolate_directions3d(cal_vals, cal_weights, cal_dirs, interp_dirs, interp_kind, smooth=0.0):
    """
    Function to parallelize interpolation. Interpolates and evaluates
    the interpolation at the new directions for one ant/pol/time etc.
    This function interpolates values on a unit sphere in 3d. This is faster because
    we can use scipy.staptial.distance.cdist euclidean for this.

    Parameters
    ----------
    cal_vals: array, values of selection, 1st axis: dir, other axes flattened
    weights: array, weights of values,
    cal_dirs: (2,n) array of floats, ra/dec values of calibrator directions
    interp_dirs: (2,n) array of floats, ra/dec values of interpolate directions
    interp_kind: string, 'wrap' or 'lin' dependent on wheter to interpolate phases
                 (independently in real and imag) or TEC, amp..
    smooth: float, smooth-parameter for interpolation. For phases, a value between 1e-3 and 0.0
        seems to work well. For TEC, less smoothing should be necessary.

    Returns
    -------
    new_vals: array of floats, interpolated values of selection at new directions
    new_weights: array of floats, weights
    """
    # approximation: Use euclidean norm in 3d instead of haversine. Every implementation of
    # haversine I found is much slower than scipy.spatial.distance.cdist for eculidean
    # to 3d for euclidean approximation. 2d-interpolating in ra-dec has the issue of wraps
    # and no cdist implementation.

    new_weights = np.ones((len(interp_dirs),cal_vals.shape[-1]))
    cal_weights[np.isnan(cal_vals)] = 0.0 # make sure all NaNs are flagged
    cal_vals[np.isnan(cal_vals)] = 0.0 # set flagged values to 0
    new_weights[:,np.sum(cal_weights, axis=0) < len(cal_weights)] # Set new weights to 0 where at least one old direction is flagged.

    if interp_kind == 'wrap': # for phases, interpolate real and imag.
        _complex = np.exp(1.j * cal_vals)
        f_interp_re = Rbf(*_cartesian_from_radec(cal_dirs), _complex.real, norm='euclidean', mode='N-D', smooth=smooth) # multiquadratic, linear
        f_interp_im = Rbf(*_cartesian_from_radec(cal_dirs), _complex.imag, norm='euclidean', mode='N-D', smooth=smooth)
        interp_re = f_interp_re(*_cartesian_from_radec(interp_dirs))
        interp_im = f_interp_im(*_cartesian_from_radec(interp_dirs))
        new_vals = np.angle(interp_re + 1.j * interp_im)
    elif interp_kind == 'lin':
        f_interp = Rbf(*_cartesian_from_radec(cal_dirs), cal_vals, norm='euclidean', mode='N-D', smooth=smooth)
        new_vals = f_interp(*_cartesian_from_radec(interp_dirs))
    return new_vals, new_weights


def run( soltab, interp_dirs, soltabOut=None, prefix='interp_', ncpu=0):
    """
    Add interpolated directions to h5parm
    Parameters
    ----------
    interp_dirs : 2d array of floats
        Shape n x 2, contains ra/dec in rad.
        For example: [[ra1,dec1],[ra2,dec2],...]

    soltabOut : string,  optional
        Default: overwrite current Soltab

    prefix : string, optional, default = "interp_".
        Name prefix of interpolated directions.
    """
    import multiprocessing as mp
    # need scipy commit f0a478c4a4172c4d2225910216de6b62721db161 for multidimensional interp. otherwise slow...
    if scipy.__version__ < '1.4.0':
        raise ImportError('SciPy version >= 1.4.0 is required to support multidimensional rbf interpolation.')

    logging.info("Working on soltab: "+soltab.name)

    # check input
    if soltab.getType() in ['phase']:
        interp_kind = 'wrap'
    elif soltab.getType() in ['tec', 'amplitude', 'rotationmeasure','tec3rd']:
        interp_kind = 'lin'
    else:
        logging.error('Soltab type {} not supported.'.format(soltab.getType()))
        return 1
    if not 'dir' in soltab.getAxesNames():
        logging.error('Data without dir axis cannot be interpolated in directions...')
        return 1
    if not 'ant' in soltab.getAxesNames():
        logging.error('Data without ant axis not supported...')
        return 1
    if len(soltab.dir) < 10:
        logging.error('Not enough directions. Use at least ten for interpolation.')
        return 1

    solset = soltab.getSolset()
    ax_ord_init = soltab.getAxesNames() # original order
    dir_ax, ant_ax = ax_ord_init.index('dir'), ax_ord_init.index('ant')
    vals, weights = soltab.val, soltab.weight
    # AXES ORDER: (dir, ant, others)
    vals, weights = np.swapaxes(vals, dir_ax, 0), np.swapaxes(weights, dir_ax, 0)
    vals, weights = np.swapaxes(vals, ant_ax, 1), np.swapaxes(weights, ant_ax, 1)
    ax_ord_swap = ax_ord_init.copy()
    ax_ord_swap[0], ax_ord_swap[dir_ax] = ax_ord_swap[dir_ax], ax_ord_swap[0]
    ax_ord_swap[1], ax_ord_swap[ant_ax] = ax_ord_swap[ant_ax], ax_ord_swap[1]
    shape_swap = vals.shape # (dir, ant, others...)

    # prepare array for interpolated values - concatenate new directions
    val_shape_inpt = list(soltab.val.shape)
    val_shape_interp = val_shape_inpt.copy()
    val_shape_interp[0] = len(interp_dirs) # match direction axis with no of interp dirs
    interp_vals = np.zeros(val_shape_interp)
    interp_weights = np.ones(val_shape_interp)


    # ra/dec of calibrator dirs. Make sure order is the same as in soltab.
    cal_dirs = np.array([soltab.getSolset().getSou()[k] for k in soltab.dir])


    # Collect arguments for parallelization on antennas
    args, selections = [], []
    # change shape: (dir, ant, others) -> flatten along others
    returnAxes = ax_ord_swap.copy()
    returnAxes.remove('ant')
    for i, (vals, weights, coord, selection) in enumerate(soltab.getValuesIter(returnAxes=returnAxes, weight=True)):
        vals = np.reshape(vals, (shape_swap[0], np.product(shape_swap[2:])))
        weights = np.reshape(weights, (shape_swap[0], np.product(shape_swap[2:])))
        args.append([vals, weights, cal_dirs, interp_dirs, interp_kind])
        selections.append(selection)

    #################### DEBUG PLOT ######################
    if False:
        # What to plot?
        antidx = 37
        sel = [500,201]
        import matplotlib.pyplot as plt
        minra, maxra = np.min(cal_dirs[:,0])-0.005, np.max(cal_dirs[:,0])+0.005
        ra_range = np.linspace(minra, maxra,500)
        mindec, maxdec = np.min(cal_dirs[:,1])-0.005, np.max(cal_dirs[:,1])+0.005
        dec_range = np.linspace(mindec, maxdec,500)
        plotdirs = np.meshgrid(ra_range,dec_range)
        plotdirs = np.array((np.array(plotdirs[0]).flatten(),np.array(plotdirs[1]).flatten())).T
        pvals, pweights = args[antidx][0:2]
        pvals, pweights = np.array(pvals)[:,sel], np.array(pweights)[:,sel]
        plotvals = np.array(interpolate_directions3d(pvals, pweights, cal_dirs, plotdirs, interp_kind, smooth=1.e-3))[0,:,0]
        plotvals = plotvals.reshape((500,500))
        fig = plt.figure(dpi=200)
        plt.xlabel('RA', labelpad=13, fontsize=7.5)
        plt.ylabel('Dec', labelpad=13, fontsize=7.5)
        cmap = plt.set_cmap('jet')
        im = plt.imshow(plotvals, cmap=cmap, extent = [minra,maxra,maxdec,mindec], vmin=-3.14, vmax=3.14)
        cb = fig.colorbar(im)
        plt.scatter(*cal_dirs.T, c=pvals[:,0], cmap=cmap, edgecolors='k', vmin=-3.14, vmax=3.14)
        plt.scatter(*cal_dirs.T, c=pvals[:,0], cmap=cmap, edgecolors='k', vmin=-3.14, vmax=3.14)
        if soltab.getType() == 'phase':
            cb.set_label('phase [rad]', fontsize=7)
        elif soltab.getType() == 'tec':
            cb.set_label('vTEC [TECU]', fontsize=7)
        cb.ax.tick_params(labelsize=6.5)
        plt.savefig('debug_screen_interp.png', dpi=200, bbox_inches='tight', )
        import sys
        sys.exit()
    # run the interpolation
    ncpu = mp.cpu_count() if ncpu == 0 else ncpu # default use all cores
    with mp.Pool(ncpu) as pool:
        logging.info('Start interpolation.')
        results = pool.starmap(interpolate_directions3d, args)

    # reorder results
    for selection, result in zip(selections,results):
        vals, weights = result
        # reconstruct shape: (dirs, ant, others)
        vals = np.reshape(vals, (val_shape_interp[dir_ax],1,*shape_swap[2:]))
        weights = np.reshape(weights, (val_shape_interp[dir_ax],1,*shape_swap[2:]))
        # reconstruct initial shape -> swap back dir and ant
        vals, weights = np.swapaxes(vals, 0, dir_ax), np.swapaxes(weights, 0, dir_ax)
        vals, weights = np.swapaxes(vals, 1, ant_ax), np.swapaxes(weights, 1, ant_ax)
        # fill output arrays
        interp_vals[tuple(selection)] = vals
        interp_weights[tuple(selection)] = weights

    # concatenate existing values and interpolated direction values
    vals = np.concatenate([soltab.val,interp_vals], axis=dir_ax)
    weights = np.concatenate([soltab.weight,interp_weights], axis=dir_ax)

    # prepare output - check if soltabOut exists
    if soltabOut in solset.getSoltabNames():
        logging.warning('Soltab {} exists. Overwriting...'.format(soltabOut))
        solset.getSoltab(soltabOut).delete()

    # set names for the interpolated directions
    interp_dir_names = np.arange(len(interp_dirs)).astype(str)
    interp_dir_names = [prefix + n.zfill(3) for n in interp_dir_names]
    # prepare axes values, append directions
    axes_vals = [soltab.getAxisValues(axisName) for axisName in ax_ord_init]
    axes_vals[dir_ax] = np.concatenate([axes_vals[dir_ax],interp_dir_names])

    # make soltabOut
    soltabout = solset.makeSoltab(soltype=soltab.getType(), soltabName=soltabOut, axesNames=ax_ord_init, \
                                  axesVals=axes_vals, vals=vals, weights=weights)

    # append interpolated dirs to solset source table
    sourceTable = solset.obj._f_get_child('source')
    # vals = [[ra, dec] for ra, dec in zip(ras, decs)]
    sourceTable.append(list(zip(*(interp_dir_names, interp_dirs))))

    # Add CREATE entry to history
    soltabout.addHistory('Created by INTERPOLATEDIRECTIONS operation from %s.' % soltab.name)
    return 0
