#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implements basic plotting

import logging
from operations_lib import *

logging.debug('Loading PLOT module.')


def make_tec_screen_plots(pp, rr, tec_fit_white, tec_fit, station_positions,
    source_names, times, height, order, beta_val, r_0, prefix = 'frame_',
    remove_gradient=True):
    """Makes plots of TEC screens"""
    import pylab
    import numpy as np
    import os
    from operations.tecscreen import calc_piercepoint
    import progressbar

    root_dir = os.path.dirname(prefix)
    if root_dir == '':
        root_dir = './'
    prestr = os.path.basename(prefix)
    try:
        os.makedirs(root_dir)
    except OSError:
        pass

    N_stations = station_positions.shape[0]
    N_sources = len(source_names)
    N_times = len(times)

    A = pylab.concatenate([pylab.kron(np.eye(N_sources),
        np.ones((N_stations,1))), pylab.kron(np.ones((N_sources,1)),
        np.eye(N_stations))], axis=1)

    N_piercepoints = N_sources * N_stations
    P = np.eye(N_piercepoints) - np.dot(np.dot(A, pylab.pinv(np.dot(A.T, A))), A.T)

    x,y,z = station_positions[0, :]
    east = np.array([-y, x, 0])
    east = east / pylab.norm(east)

    north = np.array([ -x, -y, (x*x + y*y)/z])
    north = north / pylab.norm(north)

    up = np.array([x,y,z])
    up = up / pylab.norm(up)

    T = pylab.concatenate([east[:, pylab.newaxis], north[:, pylab.newaxis]], axis=1)

    pp1 = np.dot(pp[0, :, :], T)
    lower = np.amin(pp1, axis=0)
    upper = np.amax(pp1, axis=0)
    extent = upper - lower

    lower = lower - 0.05 * extent
    upper = upper + 0.05 * extent

    extent = upper - lower

    N = 50
    xr = np.arange(lower[0], upper[0], extent[0]/N)
    yr = np.arange(lower[1], upper[1], extent[1]/N)
    screen = np.zeros((N, N, N_times))

    fitted_tec1 = tec_fit.transpose([0, 2, 1]).reshape(N_piercepoints, N_times) + np.dot(P, rr-tec_fit.transpose([0, 2, 1]).reshape(N_piercepoints, N_times))

    logging.info('Calculating TEC screen images...')
    pbar = progressbar.ProgressBar(maxval=N_times).start()
    ipbar = 0
    for k in range(N_times):
        f = tec_fit_white[:, k, :]
        for i, x in enumerate(xr[0: N]):
            for j, y in enumerate(yr[0: N]):
                p = calc_piercepoint(np.dot(np.array([x, y]), np.array([east, north])), up, height)
                d2 = np.sum(np.square(pp[k, :, :] - p[0]), axis=1)
                c = -(d2 / ( r_0**2 ) )**( beta_val / 2.0 ) / 2.0
                screen[j, i, k] = np.dot(c, f.reshape(N_piercepoints))

        # Fit and remove a gradient
        if remove_gradient:
            xs, ys = np.indices(screen.shape[0:2])
            zs = screen[:, :, k]
            XYZ = []
            for xf, yf in zip(xs.flatten().tolist(), ys.flatten().tolist()):
                XYZ.append([xf, yf, zs[xf, yf]])
            XYZ = np.array(XYZ)
            a, b, c = fitPLaneLTSQ(XYZ)
            grad_plane = a * xs + b * ys + c
            screen[:, :, k] = screen[:, :, k] - grad_plane
            for t in range(fitted_tec1.shape[0]):
                xs_pt = np.where(np.array(xr) > pp1[t, 0])[0][0]
                ys_pt = np.where(np.array(yr) > pp1[t, 1])[0][0]
                grad_plane_pt = a * ys_pt + b * xs_pt + c
                fitted_tec1[t, k] = fitted_tec1[t, k] - grad_plane_pt
        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()
    vmin = np.min([np.amin(screen), np.amin(fitted_tec1)])
    vmax = np.max([np.amax(screen), np.amax(fitted_tec1)])

    logging.info('Plotting TEC screens...')
    fig1 = pylab.figure(figsize = (7, 7))
    pbar = progressbar.ProgressBar(maxval=N_times).start()
    ipbar = 0
    for k in range(N_times):
        pylab.clf()
        im = pylab.imshow(screen[:, :, k],
            cmap = pylab.cm.jet,
            origin = 'lower',
            interpolation = 'nearest',
            extent = (xr[0], xr[-1], yr[0], yr[-1]), vmin = vmin, vmax = vmax)

        sm = pylab.cm.ScalarMappable(cmap = pylab.cm.jet,
            norm = pylab.normalize(vmin = vmin, vmax=vmax))
        sm._A = []
        pylab.title(str(i))
        cbar = pylab.colorbar()
        cbar.set_label('TECU', rotation=270)

        x = []
        y = []
        s = []
        c = []
        for j in range(fitted_tec1.shape[0]):
            x.append(pp1[j,0])
            y.append(pp1[j,1])
            xs = np.where(np.array(xr) > pp1[j,0])[0][0]
            ys = np.where(np.array(yr) > pp1[j,1])[0][0]
            s.append(max(2400*abs(fitted_tec1[j, k] - screen[ys, xs, k]), 10))
            c.append(sm.to_rgba(fitted_tec1[j, k]))

        pylab.scatter(x, y, s=s, c=c)
        labels = source_names
        for label, xl, yl in zip(labels, x[0::N_stations], y[0::N_stations]):
            pylab.annotate(
                label,
                xy = (xl, yl), xytext = (-20, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'gray', alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        pylab.title('Time {0}'.format(k))
        pylab.xlim(xr[-1], xr[0])
        pylab.ylim(yr[0], yr[-1])
        pylab.xlabel('Projected Distance (m)')
        pylab.ylabel('Projected Distance (m)')
        pylab.savefig(root_dir+'/'+prestr+'frame%0.3i.png' % k)
        pbar.update(ipbar)
        ipbar += 1
    pbar.finish()
    pylab.close(fig1)


def fitPLaneLTSQ(XYZ):
    """Fits a plane to an XYZ point cloud

    Returns (a, b, c), where Z = aX + bY + c
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

    import matplotlib.pyplot as plt
    import numpy as np
    from h5parm import solFetcher

    solsets = getParSolsets( step, parset, H )
    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    plotType = parset.getString('.'.join(["LoSoTo.Steps", step, "PlotType"]), '' )
    axesToPlot = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), '' )
    minZ, maxZ = parset.getDoubleVector('.'.join(["LoSoTo.Steps", step, "MinMax"]), [0,0] )
    prefix = parset.getString('.'.join(["LoSoTo.Steps", step, "Prefix"]), '' )

    if plotType.lower() != 'tecscreen':
        for soltab in openSoltabs( H, soltabs ):

            sf = solFetcher(soltab)
            logging.info("Plotting soltab: "+soltab._v_name)

            sf.setSelection(ant=ants, pol=pols, dir=dirs)

            # some checks
            for axis in axesToPlot:
                if axis not in sf.getAxesNames():
                    logging.error('Axis \"'+axis+'\" not found.')
                    return 1

            if (len(axesToPlot) != 2 and plotType == '2D') or \
               (len(axesToPlot) != 1 and plotType == '1D'):
                logging.error('Wrong number of axes.')
                return 1

            for vals, coord in sf.getValuesIter(returnAxes=axesToPlot):
                # TODO: implement flag control, using different color?

                title = ''
                for axis in coord:
                    if axis in axesToPlot: continue
                    title += str(coord[axis])+'_'
                title = title[:-1]

                if plotType == '2D':
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    plt.title(title)
                    plt.ylabel(axesToPlot[0])
                    plt.xlabel(axesToPlot[1])
                    p = ax.imshow(coord[axesToPlot[1]], coord[axesToPlot[0]], vals)
                    if not (minZ == 0 and maxZ == 0):
                        plt.zlim(zmin=minZ, zmax=maxZ)
                    plt.savefig(title+'.png')
                    logging.info("Saving "+prefix+title+'.png')

                if plotType == '1D':
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    plt.title(title)
                    plt.ylabel(sf.getType())
                    if not (minZ == 0 and maxZ == 0):
                        plt.ylim(ymin=minZ, ymax=maxZ)
                    plt.xlabel(axesToPlot[0])
                    p = ax.plot(coord[axesToPlot[0]], vals)
                    plt.savefig(prefix+title+'.png')
                    logging.info("Saving "+prefix+title+'.png')
    else:
        # Plot TEC screens
        i = 0
        st_tec = None
        st_pp = None
        st_tfw = None

        # Find required soltabs
        for st in openSoltabs(H, soltabs):
            if st._v_title == 'tec':
                st_tec = st
                tec_indx = i
            elif st._v_title == 'tecfitwhite':
                st_tfw = st
            elif st._v_title == 'piercepoint':
                st_pp = st
            i += 1

        if st_tec is None or st_pp is None or st_tfw is None:
            logging.warning('One or more of the required TEC solution tables '
                'not found')
            return 1

        solset = soltabs[tec_indx].split('/')[0]
        station_dict = H.getAnt(solset)
        station_names = station_dict.keys()
        station_positions = station_dict.values()
        source_dict = H.getSou(solset)
        source_names = source_dict.keys()
        source_positions = source_dict.values()

        sf_tec = solFetcher(st_tec)
        r, axis_vals = sf_tec.getValues()
        times = axis_vals['time']
        N_sources = len(source_names)
        N_times = len(times)
        N_stations = len(station_names)
        N_piercepoints = N_sources * N_stations
        rr = np.reshape(r.transpose([0, 2, 1]), [ N_piercepoints, N_times])

        sf_pp = solFetcher(st_pp)
        pp = sf_pp.getValues(retAxesVals=False)
        height = st_tfw._v_attrs['height']
        order = st_tfw._v_attrs['order']
        beta_val = st_tfw._v_attrs['beta']
        r_0 = st_tfw._v_attrs['r_0']

        sf_tfw = solFetcher(st_tfw)
        tec_fit_white = sf_tfw.getValues(retAxesVals=False)
        tec_fit = sf_tfw.getValues(weight=True, retAxesVals=False)

        make_tec_screen_plots(pp, rr, tec_fit_white, tec_fit,
            np.array(station_positions), np.array(source_names), times,
            height, order, beta_val, r_0, prefix=prefix, remove_gradient=True)

    return 0
