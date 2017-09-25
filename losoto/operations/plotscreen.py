#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation implements screen plotting
# WEIGHT: not ready

import logging
from losoto.operations_lib import *

logging.debug('Loading PLOTSCREEN module.')

def run_parser(soltab, parser, step):

    minZ, maxZ = parser.getarray( step, "MinMax", [0.0, 0.0] )
    prefix = parser.getstr( step, "Prefix", '' )
    remove_gradient = parser.getbool( step, "RemoveGradient", False )
    return run(soltab, minZ, maxZ, prefix, remove_gradient)


def normalize(phase):
    """
    Normalize phase to the range [-pi, pi].
    """
    import numpy as np

    # Convert to range [-2*pi, 2*pi].
    out = np.fmod(phase, 2.0 * np.pi)

    # Convert to range [-pi, pi]
    out[out < -np.pi] += 2.0 * np.pi
    out[out > np.pi] -= 2.0 * np.pi

    return out


def calculate_screen(inscreen, residuals, pp, N_piercepoints, k, east, north, up,
    T, Nx, Ny, sindx, height, beta_val, r_0, is_phase, outQueue):
    """
    Calculates screen images

    Parameters
    ----------
    inscreen: array
        Array of screen values at the piercepoints
    residuals: array
        Array of screen residuals at the piercepoints
    pp: array
        Array of piercepoint locations
    N_piercepoints: int
        Number of pierce points
    k: int
        Time index
    east: array
        East array
    north: array
        North array
    up: array
        Up array
    T: array
        T array
    Nx: int
        Number of pixels in x for screen
    Ny: int
        Number of pixels in y for screen
    sindx: int
        Station index
    height: float
        height of screen (m)
    beta_val: float
        power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    r_0: float
        scale size of phase fluctuations
    is_phase: bool
        input screen is a phase screen

    """
    from numpy import kron, concatenate, newaxis
    from numpy.linalg import pinv, norm
    import numpy as np
    from losoto.operations.tecscreen import calc_piercepoint

    screen = np.zeros((Nx, Ny))

    if height == 0.0:
        pp1 = pp[:, :]
    else:
        pp1 = np.dot(pp[:, :], T)

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

    D = np.resize(pp, (N_piercepoints, N_piercepoints, 3))
    D = np.transpose(D, (1, 0, 2)) - D
    D2 = np.sum(D**2, axis=2)
    C = -(D2 / r_0**2)**(beta_val / 2.0) / 2.0
    f = inscreen.reshape(N_piercepoints)
    fitted_tec = np.dot(C, f) + residuals
    if is_phase:
        fitted_tec = normalize(fitted_tec)
    for i, xi in enumerate(xr[0: Nx]):
        for j, yi in enumerate(yr[0: Ny]):
            if height == 0.0:
                p = np.array([xi, yi, 0.0])
            else:
                p, airmass = calc_piercepoint(np.dot(np.array([xi, yi]), np.array([east, north])), up, height)
            d2 = np.sum(np.square(pp - p), axis=1)
            c = -(d2 / ( r_0**2 ))**(beta_val / 2.0) / 2.0
            screen[i, j] = np.dot(c, f)

    # Calculate the piercepoint coords, normalized to 1
    if height == 0.0:
        x = pp1[:, 0]
        y = pp1[:, 1]
    else:
        x = (pp1[:, 0] / 1000.0 - min_xy[0] / 1000.0) / extent[0] # km
        y = (pp1[:, 1] / 1000.0 - min_xy[1] / 1000.0) / extent[1] # km

    outQueue.put([k, fitted_tec, screen, x, y])


def plot_frame(screen, fitted_phase1, residuals, weights, x, y, k, lower,
    upper, vmin, vmax, source_names, show_source_names, station_names, sindx,
    root_dir, prestr, is_image_plane,  midRA, midDec, outQueue):
    """
    Plots screen images

    Parameters
    ----------
    screen: array
        Array of screen values at the piercepoints
    fitted_phase1: array
        Array of fitted phase values
    residuals: array
        Array of screen residuals at the piercepoints
    weights: array
        Array of weights at the piercepoints
    x: array
        Array of piercepoint x locations
    y: array
        Array of piercepoint y locations
    k: int
        Time index
    lower: array
        Array of lower limits for plot
    upper: array
        Array of upper limits for plot
    vmin: float
        minimum value for plot range
    vmax: float
        maximum value for plot range
    source_names: list
        List of source (direction) names
    show_source_names: bool
        label sources on screen plots
    is_phase: bool

    """
    if not 'matplotlib' in sys.modules:
        import matplotlib as mpl
        mpl.rc('font',size =8 )
        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
        mpl.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt # after setting "Agg" to speed up
    from losoto.operations.phasescreen import xy2radec, makeWCS
    import numpy as np
    try:
        try:
            from astropy.visualization.wcsaxes import WCSAxes
            hasWCSaxes = True
        except:
            from wcsaxes import WCSAxes
            hasWCSaxes = True
    except:
        hasWCSaxes = False
    from matplotlib.colors import LinearSegmentedColormap


    fig = plt.figure(figsize=(7,7))

    # Set colormap
    cmap = LinearSegmentedColormap.from_list('mycmap', ['navy', 'blue',
        'turquoise', 'lightgreen', 'yellow', 'red', 'navy'])
    sm = plt.cm.ScalarMappable(cmap=cmap,
        norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []

    if is_image_plane and hasWCSaxes:
        wcs = makeWCS(midRA, midDec)
        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
    else:
        plt.gca().set_aspect('equal')
        ax = plt.gca()

    s = []
    c = []
    for j in range(fitted_phase1.shape[0]):
        if weights[j] > 0.0:
            s.append(max(10, 200*np.sqrt(weights[j]/np.median(weights))))
        else:
            s.append(10)
        c.append(sm.to_rgba(fitted_phase1[j]))

    if is_image_plane:
        min_x = np.min(x)
        max_x = np.max(x)
        min_y = np.min(y)
        max_y = np.max(y)
        extent_x = max_x - min_x
        extent_y = max_y - min_y
        lower = [min_x - 0.1 * extent_x, min_y - 0.1 * extent_y]
        upper = [max_x + 0.1 * extent_x, max_y + 0.1 * extent_y]
    else:
        # convert from m to km
        lower /= 1000.0
        upper /= 1000.0

    im = ax.imshow(screen.transpose([1, 0])[:, :],
        cmap = cmap,
        origin = 'lower',
        interpolation = 'nearest',
        extent = (lower[0], upper[0], lower[1], upper[1]),
        vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im)
    cbar.set_label('Value', rotation=270)

    ax.scatter(x, y, s=s, c=c)
    if show_source_names:
        labels = source_names
        for label, xl, yl in zip(labels, x[0::N_stations], y[0::N_stations]):
            plt.annotate(
                label,
                xy = (xl, yl), xytext = (-2, 2),
                textcoords = 'offset points', ha = 'right', va = 'bottom')

    plt.title('Station {0} ({1}), Time {2}'.format(sindx, station_names[sindx], k))
    if is_image_plane:
        ax.set_xlim(lower[0], upper[0])
        ax.set_ylim(lower[1], upper[1])
        ax.set_aspect('equal')
        if hasWCSaxes:
            RAAxis = ax.coords['ra']
            RAAxis.set_axislabel('RA', minpad=0.75)
            RAAxis.set_major_formatter('hh:mm:ss')
            DecAxis = ax.coords['dec']
            DecAxis.set_axislabel('Dec', minpad=0.75)
            DecAxis.set_major_formatter('dd:mm:ss')
            ax.coords.grid(color='black', alpha=0.5, linestyle='solid')
            plt.xlabel("RA")
            plt.ylabel("Dec")
        else:
            plt.xlabel("RA (arb. units)")
            plt.ylabel("Dec (arb. units)")
    else:
        # Reverse the axis so that RA coord increases to left
        plt.xlim(upper[0], lower[0])
        plt.ylim(lower[1], upper[1])
        plt.xlabel('Projected Distance East-West (km)')
        plt.ylabel('Projected Distance North-South (km)')
    plt.savefig(root_dir + '/' + prestr + '_station%0.4i' % sindx + '_frame%0.4i.png' % k)
    plt.close(fig)


def make_screen_plots(pp, inscreen, inresiduals, weights, station_names,
    station_positions, source_names, times, height, order, beta_val,
    r_0, prefix='frame_', remove_gradient=True, show_source_names=False,
    min_val=None, max_val=None, is_phase=False, midRA=0.0, midDec=0.0, ncpu=0):
    """
    Makes plots of screens

    Parameters
    ----------
    pp: array
        Array of piercepoint locations
    inscreen: array
        Array of screen values at the piercepoints
    residuals: array
        Array of screen residuals at the piercepoints
    weights: array
        Array of weights for each piercepoint
    source_names: array
        Array of source names
    times: array
        Array of times
    height: float
        Height of screen (m)
    order: int
        Order of screen (e.g., number of KL base vectors to keep)
    r_0: float
        Scale size of phase fluctuations
    beta_val: float
        Power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    prefix: str
        Prefix for output file names
    remove_gradient: bool
        Fit and remove a gradient from each screen
    show_source_names: bool
        Label sources on screen plots
    min_val: float
        Minimum value for plot range
    max_val: float
        Maximum value for plot range
    is_phase: bool
        Input screen is a phase screen
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees
    ncpu: int
        Number of CPUs to use

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

    # input check
    root_dir = os.path.dirname(prefix)
    if root_dir == '':
        root_dir = './'
    prestr = os.path.basename(prefix) + 'screen'
    try:
        os.makedirs(root_dir)
    except OSError:
        pass
    if ncpu == 0:
        import multiprocessing
        ncpu = multiprocessing.cpu_count()

    for sindx in range(station_positions.shape[0]):
        N_stations = 1 # screens are single-station screens
        N_sources = len(source_names)
        N_times = len(times)
        N_piercepoints = N_sources * N_stations

        xp, yp, zp = station_positions[0, :] # use first station
        east = np.array([-yp, xp, 0])
        east = east / norm(east)

        north = np.array([-xp, -yp, (xp*xp + yp*yp)/zp])
        north = north / norm(north)

        up = np.array([xp, yp, zp])
        up = up / norm(up)

        T = concatenate([east[:, newaxis], north[:, newaxis]], axis=1)

        # Use pierce point locations of first and last time slots to estimate
        # required size of plot in meters
        if height == 0.0:
            is_image_plane = True # pierce points are image plane coords
            pp1_0 = pp[0, :, 0:2]
            pp1_1 = pp[0, :, 0:2]
        else:
            is_image_plane = False
            pp1_0 = np.dot(pp[0, :, :], T)
            pp1_1 = np.dot(pp[-1, :, :], T)

        max_xy = np.amax(pp1_0, axis=0) - np.amin(pp1_0, axis=0)
        max_xy_1 = np.amax(pp1_1, axis=0) - np.amin(pp1_1, axis=0)
        if max_xy_1[0] > max_xy[0]:
            max_xy[0] = max_xy_1[0]
        if max_xy_1[1] > max_xy[1]:
            max_xy[1] = max_xy_1[1]

        min_xy = np.array([0.0, 0.0])
        extent = max_xy - min_xy
        lower = min_xy - 0.1 * extent
        upper = max_xy + 0.1 * extent
        im_extent_m = upper - lower

        residuals = inresiduals[:, :, sindx, newaxis].transpose([0, 2, 1]).reshape(N_piercepoints, N_times)
        fitted_phase1 = np.zeros((N_piercepoints, N_times))

        Nx = 24
        Ny = 0
        while Ny < 20:
            pix_per_m = Nx / im_extent_m[0]
            m_per_pix = 1.0 / pix_per_m
            Ny = int(im_extent_m[1] * pix_per_m)
            Nx += 1

        x = np.zeros((N_times, N_piercepoints)) # plot x pos of piercepoints
        y = np.zeros((N_times, N_piercepoints)) # plot y pos of piercepoints
        screen = np.zeros((Nx, Ny, N_times))

        logging.info('Calculating screen images...')
        mpm = multiprocManager(ncpu, calculate_screen)
        for k in range(N_times):
            mpm.put([inscreen[:, k, sindx], residuals[:, k], pp[k, :, :],
                N_piercepoints, k, east, north, up, T, Nx, Ny, sindx, height,
                beta_val, r_0, is_phase])
        mpm.wait()
        for (k, ft, scr, xa, ya) in mpm.get():
            screen[:, :, k] = scr
            fitted_phase1[:, k] = ft
            if is_image_plane:
                x[k, :] = xa
                y[k, :] = ya
            else:
                x[k, :] = xa - np.amin(xa) # remove offsets for each time slot
                y[k, :] = ya - np.amin(ya)

        # Normalize piercepoint locations to extent calculated above
        if not is_image_plane:
            x *= extent[0]
            y *= extent[1]

        if min_val is None:
            vmin = np.min([np.amin(screen), np.amin(fitted_phase1)])
        else:
            vmin = min_val
        if max_val is None:
            vmax = np.max([np.amax(screen), np.amax(fitted_phase1)])
        else:
            vmax = max_val

        logging.info('Plotting screens...')
        mpm = multiprocManager(ncpu, plot_frame)
        for k in range(N_times):
            mpm.put([screen[:, :, k], fitted_phase1[:, k], residuals[:, k],
            weights[:, k, sindx], x[k, :], y[k, :], k, lower, upper, vmin, vmax,
            source_names, show_source_names, station_names, sindx, root_dir,
            prestr, is_image_plane, midRA, midDec])
        mpm.wait()


def fitPLaneLTSQ(XYZ):
    """
    Fits a plane to an XYZ point cloud

    Returns (a, b, c), where Z = aX + bY + c

    Parameters
    ----------
    XYZ: array
        point cloud

    Returns
    -------
    (a, b, c): floats
        plane parameters, where Z = aX + bY + c

    """
    import numpy as np
    [rows, cols] = XYZ.shape
    G = np.ones((rows, 3))
    G[:, 0] = XYZ[:, 0]  #X
    G[:, 1] = XYZ[:, 1]  #Y
    Z = XYZ[:, 2]
    (a, b, c), resid, rank, s = np.linalg.lstsq(G, Z)
    return (a, b, c)


def run(soltab, ressoltab, minZ=-3.2, maxZ=3.2, prefix='', remove_gradient=False,
    is_phase=True, ncpu=0):
    """
    Plot screens (one plot is made per time and per station)

    Parameters
    ----------
    soltab : solution table
        Soltab containing screen
    ressoltab: solution table
        Soltab containing the screen residuals
    minZ: float, optional
        Minimum value of colorbar scale
    maxZ: float, optional
        Max value of colorbar scale
    prefix: str, optional
        String to prepend to output plots
    remove_gradient: bool, optional
        If True, remove gradient before plotting
    is_phase: bool, optional
        If True, screen is a phase screen
    ncpu: int, optional
        Number of CPUs to use. If 0, all are used

    """

    import os
    import numpy as np

    logging.info('Using input solution table: {}'.format(soltab.name))

    # Get values from soltabs
    solset = soltab.getSolset()
    screen = np.array(soltab.val)
    weights = np.array(soltab.weight)
    residuals = np.array(ressoltab.val)
    times = np.array(soltab.time)

    # Collect station and source names and positions and times, making sure
    # that they are ordered correctly.
    source_names = soltab.dir[:]
    source_dict = solset.getSou()
    source_positions = []
    for source in source_names:
        source_positions.append(source_dict[source])
    station_names = soltab.ant
    station_dict = solset.getAnt()
    station_positions = []
    for station in station_names:
        station_positions.append(station_dict[station])
    height = soltab.obj._v_attrs['height']
    order = soltab.obj._v_attrs['order']
    beta_val = soltab.obj._v_attrs['beta']
    r_0 = soltab.obj._v_attrs['r_0']
    pp = soltab.obj.piercepoint
    if height == 0.0:
        midRA = soltab.obj._v_attrs['midra']
        midDec = soltab.obj._v_attrs['middec']
    else:
        midRA = 0.0
        midDec = 0.0

    if (minZ == 0 and maxZ == 0):
        min_val = None
        max_val = None
    else:
        min_val = minZ
        max_val = maxZ

    make_screen_plots(pp, screen, residuals, weights, np.array(station_names),
        np.array(station_positions), np.array(source_names), times,
        height, order, beta_val, r_0, prefix=prefix,
        remove_gradient=remove_gradient, show_source_names=False, min_val=min_val,
        max_val=max_val, is_phase=is_phase, midRA=midRA, midDec=midDec, ncpu=ncpu)

    return 0
