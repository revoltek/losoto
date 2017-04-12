#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implements TEC screen plotting
# WEIGHT: not ready

import logging
from losoto.operations_lib import *

logging.debug('Loading PLOTTECSCREEN module.')


def make_tec_screen_plots(pp, tec_screen, residuals, station_positions,
    source_names, times, height, order, beta_val, r_0, prefix = 'frame_',
    remove_gradient=True, show_source_names=False, min_tec=None, max_tec=None):
    """Makes plots of TEC screens

    Keyword arguments:
    pp -- array of piercepoint locations
    tec_screen -- array of TEC screen values at the piercepoints
    residuals -- array of TEC screen residuals at the piercepoints
    source_names -- array of source names
    times -- array of times
    height -- height of screen (m)
    order -- order of screen (e.g., number of KL base vectors to keep)
    r_0 -- scale size of phase fluctuations (m)
    beta_val -- power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    prefix -- prefix for output file names
    remove_gradient -- fit and remove a gradient from each screen
    show_source_names -- label sources on screen plots
    min_tec -- minimum TEC value for plot range
    max_tec -- maximum TEC value for plot range
    """
    from numpy import kron, concatenate, newaxis
    from numpy.linalg import pinv, norm
    import numpy as np
    import os
    from losoto.operations.tecscreen import calc_piercepoint
    # avoids error if re-setting "agg" a second run of plot
    if not 'matplotlib' in sys.modules:
        import matplotlib as mpl
        mpl.rc('font',size =8 )
        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
        mpl.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt # after setting "Agg" to speed up
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
 
    try:
        import progressbar
    except ImportError:
        import losoto.progressbar as progressbar

    root_dir = os.path.dirname(prefix)
    if root_dir == '':
        root_dir = './'
    prestr = os.path.basename(prefix) + 'screen_'
    try:
        os.makedirs(root_dir)
    except OSError:
        pass

    N_stations = station_positions.shape[0]
    N_sources = len(source_names)
    N_times = len(times)

    A = concatenate([kron(np.eye(N_sources),
        np.ones((N_stations,1))), kron(np.ones((N_sources,1)),
        np.eye(N_stations))], axis=1)

    N_piercepoints = N_sources * N_stations
    P = np.eye(N_piercepoints) - np.dot(np.dot(A, pinv(np.dot(A.T, A))), A.T)

    x, y, z = station_positions[0, :]
    east = np.array([-y, x, 0])
    east = east / norm(east)

    north = np.array([ -x, -y, (x*x + y*y)/z])
    north = north / norm(north)

    up = np.array([x ,y, z])
    up = up / norm(up)

    T = concatenate([east[:, newaxis], north[:, newaxis]], axis=1)

    # Use pierce point locations of first time slot to estimate
    # required size of plot in meters
    pp1_0 = np.dot(pp[0, :, :], T)
    min_xy = np.amin(pp1_0, axis=0)
    max_xy = np.amax(pp1_0, axis=0)
    extent = max_xy - min_xy
    lower = min_xy - 0.1 * extent
    upper = max_xy + 0.1 * extent
    im_extent_m = upper - lower

    residuals = residuals.transpose([0, 2, 1]).reshape(N_piercepoints, N_times)
    fitted_tec1 = tec_screen.transpose([0, 2, 1]).reshape(N_piercepoints, N_times) + residuals

    Nx = 24
    Ny = 0
    while Ny < 20:
        pix_per_m = Nx / im_extent_m[0]
        m_per_pix = 1.0 / pix_per_m
        Ny = int(im_extent_m[1] * pix_per_m)
        Nx += 1

    x = [] # x coordinates for plot (km)
    y = [] # y coordinates for plot (km)
    for j in range(fitted_tec1.shape[0]):
        x.append(pp1_0[j, 0] / 1000.0)
        y.append(pp1_0[j, 1] / 1000.0)

    screen = np.zeros((Nx, Ny, N_times))
    gradient = np.zeros((Nx, Ny, N_times))

    logging.info('Calculating TEC screen images...')
    pbar = progressbar.ProgressBar(maxval=N_times).start()
    ipbar = 0
    for k in range(N_times):
        pp1 = np.dot(pp[k, :, :], T)
        min_xy = np.amin(pp1, axis=0)
        max_xy = np.amax(pp1, axis=0)
        extent = max_xy - min_xy
        lowerk = min_xy - 0.1 * extent
        upperk = max_xy + 0.1 * extent
        im_extent_mk = upperk - lowerk
        pix_per_mk = Nx / im_extent_mk[0]
        m_per_pixk = 1.0 / pix_per_mk

        xr = np.arange(lowerk[0], upperk[0], m_per_pixk)
        yr = np.arange(lowerk[1], upperk[1], m_per_pixk)

        D = np.resize(pp[k, :, :], (N_piercepoints, N_piercepoints, 3))
        D = np.transpose(D, (1, 0, 2)) - D
        D2 = np.sum(D**2, axis=2)
        C = -(D2 / r_0**2)**(beta_val / 2.0) / 2.0
        f = np.dot(pinv(C), tec_screen[:, k, :].reshape(N_piercepoints))
        for i, xi in enumerate(xr[0: Nx]):
            for j, yi in enumerate(yr[0: Ny]):
                p, airmass = calc_piercepoint(np.dot(np.array([xi, yi]), np.array([east, north])), up, height)
                d2 = np.sum(np.square(pp[k, :, :] - p), axis=1)
                c = -(d2 / ( r_0**2 ))**(beta_val / 2.0) / 2.0
                screen[i, j, k] = airmass * np.dot(c, f)

        # Fit and remove a gradient.
        if remove_gradient:
            xscr, yscr = np.indices(screen.shape[0:2])
            zs = screen[:, :, k]
            XYZ = []
            for xf, yf in zip(xscr.flatten().tolist(), yscr.flatten().tolist()):
                XYZ.append([xf, yf, zs[xf, yf]])
            XYZ = np.array(XYZ)
            a, b, c = fitPLaneLTSQ(XYZ)
            grad_plane = a * xscr + b * yscr + c
            gradient[:, :, k] = grad_plane
            screen[:, :, k] = screen[:, :, k] - grad_plane
            screen[:, :, k] = screen[:, :, k] - np.mean(screen[:, :, k])

            # Match fitted values to gradient-free screen
            for t in range(fitted_tec1.shape[0]):
                xs_pt = (x[t] * 1000.0 - lower[0]) / m_per_pix
                ys_pt = (y[t] * 1000.0 - lower[1]) / m_per_pix
                fitted_tec1[t, k] = screen[xs_pt, ys_pt, k] + residuals[t, k]

        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()
    if min_tec is None:
        vmin = np.min([np.amin(screen), np.amin(fitted_tec1)])
    else:
        vmin = min_tec
    if max_tec is None:
        vmax = np.max([np.amax(screen), np.amax(fitted_tec1)])
    else:
        vmax = max_tec

    logging.info('Plotting TEC screens...')
    fig, ax = plt.subplots(figsize=[7, 7])
    pbar = progressbar.ProgressBar(maxval=N_times).start()
    ipbar = 0
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet,
        norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    plt.gca().set_aspect('equal')

    for k in range(N_times):
        s = []
        c = []
        for j in range(fitted_tec1.shape[0]):
            fit_screen_diff = abs(residuals[j, k])
            s.append(max(20*fit_screen_diff/0.01, 10))
            c.append(sm.to_rgba(fitted_tec1[j, k]))

        plt.clf()
        im = plt.imshow(screen.transpose([1, 0, 2])[:, :, k],
            cmap = plt.cm.jet,
            origin = 'lower',
            interpolation = 'nearest',
            extent = (lower[0]/1000.0, upper[0]/1000.0, lower[1]/1000.0, upper[1]/1000.0),
            vmin=vmin, vmax=vmax)

        cbar = plt.colorbar(im)
        cbar.set_label('TECU', rotation=270)

        plt.scatter(x, y, s=s, c=c)
        if show_source_names:
            labels = source_names
            for label, xl, yl in zip(labels, x[0::N_stations], y[0::N_stations]):
                plt.annotate(
                    label,
                    xy = (xl, yl), xytext = (-2, 2),
                    textcoords = 'offset points', ha = 'right', va = 'bottom')

        plt.title('Screen {0}'.format(k))
        plt.xlim(upper[0]/1000.0, lower[0]/1000.0)
        plt.ylim(lower[1]/1000.0, upper[1]/1000.0)
        plt.xlabel('Projected Distance along RA (km)')
        plt.ylabel('Projected Distance along Dec (km)')

        if remove_gradient:
            pp1 = np.dot(pp[k, :, :], T)
            min_xy = np.amin(pp1, axis=0)
            max_xy = np.amax(pp1, axis=0)
            extent = max_xy - min_xy
            lowerk = min_xy - 0.05 * extent
            upperk = max_xy + 0.05 * extent
            xr = np.arange(lowerk[0], upperk[0], m_per_pix)
            yr = np.arange(lowerk[1], upperk[1], m_per_pix)
            axins = inset_axes(ax, width="15%", height="10%", loc=2)
            axins.imshow(gradient.transpose([1, 0, 2])[:, : ,k],
                cmap = plt.cm.jet,
                origin = 'lower',
                interpolation = 'nearest',
                extent = (xr[0]/1000.0, xr[-1]/1000.0, yr[0]/1000.0, yr[-1]/1000.0),
                vmin=vmin, vmax=vmax)
            plt.xticks(visible=False)
            plt.yticks(visible=False)
            axins.set_xlim(xr[-1]/1000.0, xr[0]/1000.0)
            axins.set_ylim(yr[0]/1000.0, yr[-1]/1000.0)

        plt.savefig(root_dir+'/'+prestr+'frame%0.4i.png' % k)
        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()
    plt.close(fig)


def fitPLaneLTSQ(XYZ):
    """Fits a plane to an XYZ point cloud

    Returns (a, b, c), where Z = aX + bY + c

    Keyword arguments:
    XYZ -- point cloud
    """
    import numpy as np
    [rows, cols] = XYZ.shape
    G = np.ones((rows, 3))
    G[:, 0] = XYZ[:, 0]  #X
    G[:, 1] = XYZ[:, 1]  #Y
    Z = XYZ[:, 2]
    (a, b, c), resid, rank, s = np.linalg.lstsq(G, Z)
    return (a, b, c)


def run( step, parset, H ):

    import os
    import numpy as np
    from losoto.h5parm import solFetcher, solHandler

    soltabs = getParSoltabs( step, parset, H )

    minZ, maxZ = parset.getDoubleVector('.'.join(["LoSoTo.Steps", step, "MinMax"]), [0,0] )
    prefix = parset.getString('.'.join(["LoSoTo.Steps", step, "Prefix"]), '' )
    remove_gradient = parset.getBool('.'.join(["LoSoTo.Steps", step, "RemoveGradient"]), False )

    # Plot various TEC-screen properties

    for st_scr in openSoltabs(H, soltabs):

        # Check if soltab is a tecscreen table
        full_name = st_scr._v_parent._v_name + '/' + st_scr._v_name
        if st_scr._v_title != 'tecscreen':
            logging.warning('Solution table {0} is not a tecscreen solution '
                'table. Skipping.'.format(full_name))
            continue
        logging.info('Using input solution table: {0}'.format(full_name))

        # Plot TEC screens as images
        solset = st_scr._v_parent
        sf_scr = solFetcher(st_scr)
        r, axis_vals = sf_scr.getValues()
        source_names = axis_vals['dir']
        station_names = axis_vals['ant']
        station_dict = H.getAnt(solset)
        station_positions = []
        for station in station_names:
            station_positions.append(station_dict[station])
        times = axis_vals['time']

        tec_screen, axis_vals = sf_scr.getValues()
        times = axis_vals['time']
        residuals = sf_scr.getValues(weight=True, retAxesVals=False)
        height = st_scr._v_attrs['height']
        order = st_scr._v_attrs['order']
        beta_val = st_scr._v_attrs['beta']
        r_0 = st_scr._v_attrs['r_0']
        pp = sf_scr.t.piercepoint

        if (minZ == 0 and maxZ == 0):
            min_tec = None
            max_tec = None
        else:
            min_tec = minZ
            max_tec = maxZ

        make_tec_screen_plots(pp, tec_screen, residuals,
            np.array(station_positions), np.array(source_names), times,
            height, order, beta_val, r_0, prefix=prefix,
            remove_gradient=remove_gradient, show_source_names=False, min_tec=min_tec,
            max_tec=max_tec)

    return 0
