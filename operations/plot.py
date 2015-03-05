#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implements basic plotting
# WEIGHT: flag-only compliant

import logging
from operations_lib import *

logging.debug('Loading PLOT module.')

def run( step, parset, H ):

    import os
    import numpy as np
    from itertools import cycle
    from h5parm import solFetcher, solHandler
    # avoids error if re-setting "agg" a second run of plot
    if not 'matplotlib' in sys.modules:
        import matplotlib as mpl
        mpl.rc('font',size =8 )
        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
        mpl.use("Agg")
    import matplotlib.pyplot as plt # after setting "Agg" to speed up

    soltabs = getParSoltabs( step, parset, H )

    # 1D or 2D array in form X, [Y]
    axesInPlot = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    minZ, maxZ = parset.getDoubleVector('.'.join(["LoSoTo.Steps", step, "MinMax"]), [0,0] )
    # the axis to plot on one page - e.g. ant to get all antenna's on one plot #
    axisInTable = parset.getString('.'.join(["LoSoTo.Steps", step, "TableAxis"]), '' )
    # the axis to plot in different colours - e.g. pol to get all correlations on one plot #
    axisInCol = parset.getString('.'.join(["LoSoTo.Steps", step, "ColorAxis"]), '' )
    # the axis to plot is different shades (alpha) - e.g. freq for a small range to compare subband to subband solutions on one plot #
    axisInShade = parset.getString('.'.join(["LoSoTo.Steps", step, "ShadeAxis"]), '' )
    # log='XYZ' to set which axes to put in Log
    log = parset.getString('.'.join(["LoSoTo.Steps", step, "Log"]), "" )
    plotflagged = parset.getBool('.'.join(["LoSoTo.Steps", step, "PlotFLagged"]), False )
    dounwrap = parset.getBool('.'.join(["LoSoTo.Steps", step, "Unwrap"]), False )
    prefix = parset.getString('.'.join(["LoSoTo.Steps", step, "Prefix"]), '' )

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Plotting soltab: "+soltab._v_name)
        sf = solFetcher(soltab)

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        # some checks
        for axis in axesInPlot:
            if axis not in sf.getAxesNames():
                logging.error('Axis \"'+axis+'\" not found.')
                return 1

        if len(axisInPlot) ==2 :
            3D = True
            # not color/shade possible in 3D
            axisInCol == []
            axisInShade == []
        elif len(axisInPlot) != 1:
            3D = False
            logging.error('Axes must be a len 1 or 2 array.')
            return 1

        if len(set(axisInTable+axesInPlot+axisInCol+axisInShade)) != len(axisInTable+axesInPlot+axisInCol+axisInShade):
            logging.error('Axis defined multiple times.')
            return 1

        # all axes that are not iterated by anything else
        axesInFile = sf.getAxesNames()
        for axis in axisInTable+axesInPlot+axisInCol+axisInShade:
            axesInFile.remove(axis)
 
        # set subplots scheme
        if axisInTable != '':
            Nplots = sf.getAxisLen(axisInTable)
        else:
            Nplots = 1
        Nr = int(np.ceil(np.sqrt(Nplots)))
        Nc = int(np.ceil(np.float(Nplots)/Nr))

        # cycle on files
        for vals, coord, selection in sf.getValuesIter(returnAxes=axisInTable+axisInCol+axisInShade+axesInPlot):
            
            fig = plt.figure()

            # set filename
            filename = ''
            for axis in axesInFile:
                filename += axis+''+str(coord[axis])+'_'
            filename = filename[:-1] # remove last _

            # create multiplot
            figgrid, axa = plt.subplots(Nr, Nc, figsize=(16,12), sharex=True, sharey=True)

            sf.selection = selection
            # cycle on tables
            for vals, coord, selection in sf.getValuesIter(returnAxes=axisInCol+axisInShade+axesInPlot):

                # set log scales if activated
                if 'X' in log: axa.set_xscale('log')
                if 'Y' in log: axa.set_yscale('log')
                if 3D and 'Z' in log: axa.set_zscale('log')

                # set tile
                title = ''
                for axis in coord:
                    if axis in axesInPlot+axisInCol+axisInShade: continue
                    title += axis+':'+str(coord[axis])+' '
                title = title[:-1] # remove last ' '

                # axes label and values
                xvals = coord[axesInPlot[0]]
                axa[Nr-1][0].set_xlabel(axesToPlot[0])
                if 3D:
                    yvals = coord[axesInPlot[1]]
                    axa[Nr-1][0].set_ylabel(axesToPlot[1])
                    axa[Nr-1][0].set_zlabel(sf.getType())
                else:
                    axa[Nr-1][0].set_ylabel(sf.getType())
                
                # if plotting antenna - convert to number
                if axesToPlot[0] == 'ant':
                    xvals = np.arange(len(xvals))
                # if plotting time - convert in h/min/s
                elif axesToPlot[0] == 'time':
                    if len(xvals) > 3600:
                        xvals = (xvals-xvals[0])/3600.  # hrs
                        xlabelunit = '[hr]'
                    elif len(xvals) > 60:
                        xvals = (xvals-xvals[0])/60.   # mins
                        xlabelunit = '[min]'
                    else:
                        xvals = (xvals-xvals[0])  # sec
                        xlabelunit = '[s]'
                    axa[Nr-1][0].set_xlabel(axesToPlot[0]+' '+xlabelunit)
                # if plotting freq convert in MHz
                elif axesToPlot[0] == 'freq': 
                    xvals = xvals/1.e6 # MHz
                    xlabelunit = '[MHz]'
                    axa[Nr-1][0].set_xlabel(axesToPlot[0]+' '+xlabelunit)

                # same as above but for y-axis
                if 3D and xesToPlot[1] == 'ant':
                    yvals = np.arange(len(yvals))
                elif 3D and axesToPlot[1] == 'time':
                    if len(xvals) > 3600:
                        yvals = (yvals-yvals[0])/3600.  # hrs
                        ylabelunit = '[hr]'
                    elif len(xvals) > 60:
                        yvals = (yvals-yvals[0])/60.   # mins
                        ylabelunit = '[min]'
                    else:
                        yvals = (yvals-yvals[0])  # sec
                        ylabelunit = '[s]'
                    axa[Nr-1][0].set_ylabel(axesToPlot[1]+' '+ylabelunit)
                elif axesToPlot[0] == 'freq':  # Mhz
                    yvals = yvals/1.e6
                    xlabelunit = '[MHz]'
                    axa[Nr-1][0].set_ylabel(axesToPlot[1]+' '+ylabelunit)

                # set colors (red reserved for flags)
                colors = cycle(['g', 'b', 'c', 'm', 'y', 'k'])

                sf.selection = selection
                # cycle on colors
                for vals, coord, selection in sf.getValuesIter(returnAxes=axisInShade+axesInPlot):

                    # set color
                    color = next(colors)

                    # set shades
                    if axisInShade != '':
                        shades = cycle(np.arange(0.1,1,0.9/sf.getAxesLen(axisInShade)))
                    else
                        shades = cycle([1])

                    sf.selection = selection
                    # cycle on shades
                    for vals, coord, selection in sf.getValuesIter(returnAxes=axesInPlot):

                        # set shade
                        shade = next(shades)

                        sf.selection = selection
                        # finally cycle on lines
                        for vals, weight, coord, selection in sf.getValuesIter(returnAxes=axesToPlot, weight=True):

                        # unwrap if required
                        if (sf.getType() == 'phase' or sf.getType() == 'scalarphase') and dounwrap: vals = unwrap(vals)
            
                        # plotting
                        if 3D:
                            if log: plt.pcolormesh(xvals, yvals , np.log10(vals))
                            else: plt.pcolormesh(xvals, yvals, vals)
                            if not (minZ == 0 and maxZ == 0):
                                plt.zlim(zmin=minZ, zmax=maxZ)
                            plt.colorbar()
                        else:
                            if sf.getType() == 'amplitude':
                                p = ax.plot(xvals[np.where(weight!=0)], vals[np.where(weight!=0)], 'k-', color=color)
                            else:
                                p = ax.plot(xvals[np.where(weight!=0)], vals[np.where(weight!=0)], 'k.', color=color)
                            if plotflagged: p = ax.plot(xvals[np.where(weight==0)], vals[np.where(weight==0)], 'ro') # plot flagged points
                            if not (minZ == 0 and maxZ == 0):
                                plt.ylim(ymin=minZ, ymax=maxZ)

            logging.info("Saving "+prefix+title+'.png')
            try:
                plt.savefig(prefix+title+'.png')
            except:
                logging.error('Error saving file, wrong path?')
                return 1
            plt.close(fig)

    return 0
