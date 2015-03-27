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
    from itertools import cycle, chain
    from h5parm import solFetcher, solHandler
    # avoids error if re-setting "agg" a second run of plot
    if not 'matplotlib' in sys.modules:
        import matplotlib as mpl
        mpl.rc('font',size =8 )
        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
        mpl.use("Agg")
    import matplotlib.pyplot as plt # after setting "Agg" to speed up

    soltabs = getParSoltabs( step, parset, H )

    # 1- or 2-element array in form X, [Y]
    axesInPlot = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    minZ, maxZ = parset.getDoubleVector('.'.join(["LoSoTo.Steps", step, "MinMax"]), [0,0] )
    # the axis to plot on one page - e.g. ant to get all antenna's on one plot #
    axisInTable = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "TableAxis"]), [] )
    # the axis to plot in different colours - e.g. pol to get all correlations on one plot #
    axisInCol = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "ColorAxis"]), [] )
    # the axis to plot is different shades (alpha) - e.g. freq for a small range to compare subband to subband solutions on one plot #
    axisInShade = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "ShadeAxis"]), [] )
    # log='XYZ' to set which axes to put in Log
    log = parset.getString('.'.join(["LoSoTo.Steps", step, "Log"]), "" )
    plotflagged = parset.getBool('.'.join(["LoSoTo.Steps", step, "PlotFlagged"]), False )
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

        cmesh = False
        if len(axesInPlot) == 2:
            cmesh = True
            # not color/shade possible in 3D
            axisInCol == []
            axisInShade == []
        elif len(axesInPlot) != 1:
            logging.error('Axes must be a len 1 or 2 array.')
            return 1

        if len(set(axisInTable+axesInPlot+axisInCol+axisInShade)) != len(axisInTable+axesInPlot+axisInCol+axisInShade):
            logging.error('Axis defined multiple times.')
            return 1

        if len(axisInTable) > 1 or len(axisInCol) > 1 or len(axisInShade) > 1:
            logging.error('Too many AxisInTable/AxisInCol/AxisInShade, they must be at most one each.')
            return 1

        # all axes that are not iterated by anything else
        axesInFile = sf.getAxesNames()
        for axis in axisInTable+axesInPlot+axisInCol+axisInShade:
            axesInFile.remove(axis)
 
        # set subplots scheme
        if axisInTable != []:
            Nplots = sf.getAxisLen(axisInTable[0])
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
            print "Filename", filename

            # create multiplot
            figgrid, axa = plt.subplots(Nr, Nc, figsize=(16,12), sharex=True, sharey=True)
            figgrid.subplots_adjust(hspace=0, wspace=0)
            axaiter = chain.from_iterable(axa)

            # axis vals (they are always the same, regulat arrays)
            xvals = coord[axesInPlot[0]]
            # if plotting antenna - convert to number
            if axesInPlot[0] == 'ant':
                xvals = np.arange(len(xvals))
            
            # if plotting time - convert in h/min/s
            xlabelunit=''
            if axesInPlot[0] == 'time':
                if len(xvals) > 3600:
                    xvals = (xvals-xvals[0])/3600.  # hrs
                    xlabelunit = ' [hr]'
                elif len(xvals) > 60:
                    xvals = (xvals-xvals[0])/60.   # mins
                    xlabelunit = ' [min]'
                else:
                    xvals = (xvals-xvals[0])  # sec
                    xlabelunit = ' [s]'
            # if plotting freq convert in MHz
            elif axesInPlot[0] == 'freq': 
                xvals = xvals/1.e6 # MHz
                xlabelunit = ' [MHz]'

            if cmesh:
                # axis vals (they are always the same, regulat arrays)
                yvals = coord[axesInPlot[1]]
                # same as above but for y-axis
                if cmesh and axesInPlot[1] == 'ant':
                    yvals = np.arange(len(yvals))

                ylabelunit=''
                if axesInPlot[1] == 'time':
                    if len(xvals) > 3600:
                        yvals = (yvals-yvals[0])/3600.  # hrs
                        ylabelunit = ' [hr]'
                    elif len(xvals) > 60:
                        yvals = (yvals-yvals[0])/60.   # mins
                        ylabelunit = ' [min]'
                    else:
                        yvals = (yvals-yvals[0])  # sec
                        ylabelunit = ' [s]'
                elif axesInPlot[0] == 'freq':  # Mhz
                    yvals = yvals/1.e6
                    ylabelunit = ' [MHz]'

            # axes label 
            [ax.set_xlabel(axesInPlot[0]+xlabelunit) for ax in axa[-1,:]]
            if cmesh:
                [ax.set_ylabel(axesInPlot[1]+ylabelunit) for ax in axa[:,0]]
            else:
                [ax.set_ylabel(sf.getType()) for ax in axa[:,0]]

            sf2 = solFetcher(soltab)
            sf2.selection = selection
            # cycle on tables
            for vals, coord, selection in sf2.getValuesIter(returnAxes=axisInCol+axisInShade+axesInPlot):

                # this axa
                ax = next(axaiter)

                # set tile
                title = ''
                for axis in coord:
                    if axis in axesInPlot+axisInCol+axisInShade: continue
                    title += axis+':'+str(coord[axis])+' '
                title = title[:-1] # remove last ' '
                ax.text(.5, .9, title, horizontalalignment='center',fontsize=8,transform=ax.transAxes)
               
                # set log scales if activated
                if 'X' in log: ax.set_xscale('log')
                if 'Y' in log: ax.set_yscale('log')

                # set colors (red reserved for flags)
                colors = cycle(['g', 'b', 'c', 'm', 'y', 'k'])

                sf3 = solFetcher(soltab)
                sf3.selection = selection
                # cycle on colors
                for vals, coord, selection in sf3.getValuesIter(returnAxes=axisInShade+axesInPlot):

                    # set color
                    color = next(colors)

                    # set shades
                    if axisInShade != []:
                        shades = cycle(np.arange(0.1,1,0.9/sf.getAxisLen(axisInShade[0])))
                    else:
                        shades = cycle([1])

                    sf4 = solFetcher(soltab)
                    sf4.selection = selection
                    # cycle on shades
                    for vals, weight, coord, selection in sf4.getValuesIter(returnAxes=axesInPlot, weight=True):

                        # set shade
                        shade = next(shades)

                        # unwrap if required
                        if (sf.getType() == 'phase' or sf.getType() == 'scalarphase') and dounwrap:
                            vals = unwrap(vals)
        
                        # plotting
                        if cmesh:
                            if log: ax.pcolormesh(xvals, yvals , np.log10(vals))
                            else: ax.pcolormesh(xvals, yvals, vals)
                            if not (minZ == 0 and maxZ == 0):
                                ax.zlim(zmin=minZ, zmax=maxZ)
                            plt.colorbar(label=sf.getType())
                        else:
                            if sf.getType() == 'amplitude':
                                ax.plot(xvals[np.where(weight!=0)], vals[np.where(weight!=0)], 'k-', color=color)
                            else:
                                ax.plot(xvals[np.where(weight!=0)], vals[np.where(weight!=0)], 'k.', color=color)
                            if plotflagged: ax.plot(xvals[np.where(weight==0)], vals[np.where(weight==0)], 'ro') # plot flagged points
                            if not (minZ == 0 and maxZ == 0):
                                plt.ylim(ymin=minZ, ymax=maxZ)

            logging.info("Saving "+prefix+filename+'.png')
            try:
                plt.savefig(prefix+filename+'.png', bbox_inches='tight')
            except:
                logging.error('Error saving file, wrong path?')
                return 1
            plt.close(fig)

    return 0
