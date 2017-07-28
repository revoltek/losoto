#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implements basic plotting
# WEIGHT: flag-only compliant

import logging
from losoto.operations_lib import *

logging.debug('Loading PLOT module.')

def plot(Nplots, NColFig, figSize, cmesh, axesInPlot, axisInTable, xvals, yvals, xlabelunit, ylabelunit, datatype, filename, titles, log, dataCube, minZ, maxZ, plotflag, makeMovie, antCoords, outQueue):
        import os
        from itertools import cycle, chain
        import numpy as np
        # avoids error if re-setting "agg" a second run of plot
        if not 'matplotlib' in sys.modules:
            import matplotlib as mpl
            mpl.rcParams['xtick.labelsize'] = 20
            mpl.rcParams['font.size'] = 20
            #mpl.rc('figure.subplot',left=0.1, bottom=0.1, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
            mpl.use("Agg")
        import matplotlib.pyplot as plt # after setting "Agg" to speed up


        autominZ = np.inf; automaxZ = -np.inf

        # if user-defined number of col use that
        if NColFig != 0: Nc = NColFig
        else: Nc = int(np.ceil(np.sqrt(Nplots)))
        Nr = int(np.ceil(np.float(Nplots)/Nc))

        if figSize[0] == 0:
            if makeMovie: figSize[0]=5+2*Nc
            else: figSize[0]=10+3*Nc
        if figSize[1] == 0:
            if makeMovie: figSize[1]=4+1*Nr
            else: figSize[1]=8+2*Nr
        
        figgrid, axa = plt.subplots(Nr, Nc, figsize=figSize, sharex=True, sharey=True)

        if Nplots == 1: axa = np.array([axa])
        figgrid.subplots_adjust(hspace=0, wspace=0)
        axaiter = chain.from_iterable(axa)

        # axes label 
        if len(axa.shape) == 1: # only one row 
            [ax.set_xlabel(axesInPlot[0]+xlabelunit, fontsize=20) for ax in axa[:]]
            if cmesh:
                axa[0].set_ylabel(axesInPlot[1]+ylabelunit, fontsize=20)
            else:
                axa[0].set_ylabel(datatype, fontsize=20)
        else:
            [ax.set_xlabel(axesInPlot[0]+xlabelunit, fontsize=20) for ax in axa[-1,:]]
            if cmesh:
                [ax.set_ylabel(axesInPlot[1]+ylabelunit, fontsize=20) for ax in axa[:,0]]
            else:
                [ax.set_ylabel(datatype, fontsize=20) for ax in axa[:,0]]

        for Ntab, title in enumerate(titles):
           
            ax = axa.flatten()[Ntab]
            ax.text(.5, .9, title, horizontalalignment='center', fontsize=14, transform=ax.transAxes)
           
            # set log scales if activated
            if 'X' in log: ax.set_xscale('log')
            if 'Y' in log: ax.set_yscale('log')

            colors = cycle(['g', 'b', 'c', 'm'])
            for Ncol, data in enumerate(dataCube[Ntab]):

                # set color, use defined colors if a few lines, otherwise a continuum colormap
                if len(dataCube[Ntab]) <= 4:
                    color = colors.next()
                    colorFlag = 'r'
                else:
                    color = plt.cm.jet(Ncol/float(len(dataCube[Ntab])-1)) # from 0 to 1
                    colorFlag = 'k'
                vals = dataCube[Ntab][Ncol]

                # plotting
                if cmesh:
                    # setting min max
                    if minZ == None and maxZ == None: 
                        autominZ = np.mean(vals) - 3*np.std(vals)
                        automaxZ = np.mean(vals) + 3*np.std(vals)
                    else:
                        autominZ = minZ
                        automaxZ = maxZ
                   # stratch the imshow output to fill the plot size
                    bbox = ax.get_window_extent().transformed(figgrid.dpi_scale_trans.inverted())
                    aspect = ((xvals[-1]-xvals[0])*bbox.height)/((yvals[-1]-yvals[0])*bbox.width)
                    if 'Z' in log:
                        if minZ != None:
                            minZ = np.log10(minZ)
                        if maxZ != None:
                            maxZ = np.log10(maxZ)
                        vals = np.log10(vals)
                    ax.imshow(vals, origin='lower', interpolation="none", cmap=plt.cm.jet, extent=[xvals[0],xvals[-1],yvals[0],yvals[-1]], aspect=str(aspect), vmin=autominZ, vmax=automaxZ)
                # make an antenna plot
                elif antCoords != []:
                    ax.set_xlabel('')
                    ax.set_ylabel('')
                    ax.axes.get_xaxis().set_ticks([])
                    ax.axes.get_yaxis().set_ticks([])
                    areas = 15 + np.pi * (10 * ( vals+np.abs(np.min(vals)) ) / np.max( vals+np.abs(np.min(vals)) ))**2 # normalize marker diameter to 15-30 pt
                    ax.scatter(antCoords[0], antCoords[1], c=vals, s=areas)
                else:
                    ax.plot(xvals, vals, 'o', color=color, markersize=2, markeredgecolor='none') # flagged data are automatically masked
                    if plotflag: 
                        ax.plot(xvals[vals.mask], vals.data[vals.mask], 'o', color=colorFlag, markersize=2, markeredgecolor='none') # plot flagged points
                    ax.set_xlim(xmin=min(xvals), xmax=max(xvals))

                    # find proper min max as the automatic setting is shit
                    if not all(vals.mask == True):
                        if autominZ > vals.min(fill_value=np.inf) or autominZ == np.inf: 
                            autominZ = vals.min(fill_value=np.inf)
                        if automaxZ < vals.max(fill_value=-np.inf) or automaxZ == -np.inf:
                            automaxZ = vals.max(fill_value=-np.inf)

        if not cmesh: 
            if minZ is not None:
                ax.set_ylim(ymin=minZ)
            else:
                ax.set_ylim(ymin=autominZ)
            if maxZ is not None:
                ax.set_ylim(ymax=maxZ)
            else:
                ax.set_ylim(ymax=automaxZ)

        logging.info("Saving "+filename+'.png')
        try:
            figgrid.savefig(filename+'.png', bbox_inches='tight')
        except:
            figgrid.tight_layout()
            figgrid.savefig(filename+'.png')
        plt.close()


def run( step, parset, H ):

    import os, random
    import numpy as np
    from losoto.h5parm import solFetcher, solHandler

    def normalize(phase):
        """
        Normalize phase to the range [-pi, pi].
        """
        # Convert to range [-2*pi, 2*pi].
        out = np.fmod(phase, 2.0 * np.pi)
        # Remove nans
        np.putmask(out, out!=out, 0)
        # Convert to range [-pi, pi]
        out[out < -np.pi] += 2.0 * np.pi
        out[out > np.pi] -= 2.0 * np.pi
        return out

    soltabs = getParSoltabs( step, parset, H )

    # 1- or 2-element array in form X, [Y]
    axesInPlot = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    NColFig = parset.getInt('.'.join(["LoSoTo.Steps", step, "Columns"]), 0 )
    figSize = parset.getIntVector('.'.join(["LoSoTo.Steps", step, "FigSize"]), [0,0] )
    minZ, maxZ = parset.getDoubleVector('.'.join(["LoSoTo.Steps", step, "MinMax"]), [0,0] )
    if minZ == 0: minZ = None
    if maxZ == 0: maxZ = None
    axisInTable = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "TableAxis"]), [] )
    axisInCol = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "ColorAxis"]), [] )
    axisInDiff = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "DiffAxis"]), [] )
    # log='XYZ' to set which axes to put in Log
    log = parset.getString('.'.join(["LoSoTo.Steps", step, "Log"]), "" )
    plotflag = parset.getBool('.'.join(["LoSoTo.Steps", step, "PlotFlag"]), False )
    dounwrap = parset.getBool('.'.join(["LoSoTo.Steps", step, "Unwrap"]), False )
    ref = parset.getString('.'.join(["LoSoTo.Steps", step, "Reference"]), '' )
    tablesToAdd = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Add"]), [] )
    makeAntPlot = parset.getBool('.'.join(["LoSoTo.Steps", step, "MakeAntPlot"]), False )
    makeMovie = parset.getBool('.'.join(["LoSoTo.Steps", step, "MakeMovie"]), False )
    prefix = parset.getString('.'.join(["LoSoTo.Steps", step, "Prefix"]), '' )

    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 0 )
    if ncpu == 0:
        import multiprocessing
        ncpu = multiprocessing.cpu_count()

    if makeMovie: 
        prefix = prefix+'__tmp__'

    if os.path.dirname(prefix) != '' and not os.path.exists(os.path.dirname(prefix)):
        logging.debug('Creating '+os.path.dirname(prefix)+'.')
        os.makedirs(os.path.dirname(prefix))

    if ref == '': ref = None
    sfsAdd = [ solFetcher(soltab) for soltab in openSoltabs(H, tablesToAdd) ]

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
            # not color possible in 3D
            axisInCol = []
        elif len(axesInPlot) != 1:
            logging.error('Axes must be a len 1 or 2 array.')
            return 1

        if len(set(axisInTable+axesInPlot+axisInCol+axisInDiff)) != len(axisInTable+axesInPlot+axisInCol+axisInDiff):
            logging.error('Axis defined multiple times.')
            return 1

        # just because we use lists, check that they are 1-d
        if len(axisInTable) > 1 or len(axisInCol) > 1 or len(axisInDiff) > 1:
            logging.error('Too many TableAxis/ColAxis/DiffAxis, they must be at most one each.')
            return 1

        # all axes that are not iterated by anything else
        axesInFile = sf.getAxesNames()
        for axis in axisInTable+axesInPlot+axisInCol+axisInDiff:
            axesInFile.remove(axis)
 
        # set subplots scheme
        if axisInTable != []:
            Nplots = sf.getAxisLen(axisInTable[0])
        else:
            Nplots = 1

        # prepare antennas coord in makeAntPlot case
        if makeAntPlot:
            if axesInPlot != ['ant']:
                logging.error('If makeAntPlot is selected the "Axes" values must be "ant"')
                return 1
            antCoords = [[],[]]
            for ant in sf.getAxisValues('ant'): # select only user-selected antenna in proper order
                antCoords[0].append(H.getAnt(sf.getAddress().split('/')[0])[ant][0])
                antCoords[1].append(H.getAnt(sf.getAddress().split('/')[0])[ant][1])
        else:
            antCoords = []
            
        datatype = sf.getType()

        # start processes for multi-thread
        mpm = multiprocManager(ncpu, plot)

        # cycle on files
        if makeMovie: pngs = [] # store png filenames
        for vals, coord, selection in sf.getValuesIter(returnAxes=axisInDiff+axisInTable+axisInCol+axesInPlot):
           
            # set filename
            filename = ''
            for axis in axesInFile:
                filename += axis+str(coord[axis])+'_'
            filename = filename[:-1] # remove last _

            # axis vals (they are always the same, regulat arrays)
            xvals = coord[axesInPlot[0]]
            # if plotting antenna - convert to number
            if axesInPlot[0] == 'ant':
                xvals = np.arange(len(xvals))
            
            # if plotting time - convert in h/min/s
            xlabelunit=''
            if axesInPlot[0] == 'time':
                if xvals[-1] - xvals[0] > 3600:
                    xvals = (xvals-xvals[0])/3600.  # hrs
                    xlabelunit = ' [hr]'
                elif xvals[-1] - xvals[0] > 60:
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
                # axis vals (they are always the same, regular arrays)
                yvals = coord[axesInPlot[1]]
                # same as above but for y-axis
                if axesInPlot[1] == 'ant':
                    yvals = np.arange(len(yvals))

                if len(xvals) <= 1 or len(yvals) <=1:
                    logging.error('3D plot must have more then one value per axes.')
                    mpm.wait()
                    return 1

                ylabelunit=''
                if axesInPlot[1] == 'time':
                    if yvals[-1] - yvals[0] > 3600:
                        yvals = (yvals-yvals[0])/3600.  # hrs
                        ylabelunit = ' [hr]'
                    elif yvals[-1] - yvals[0] > 60:
                        yvals = (yvals-yvals[0])/60.   # mins
                        ylabelunit = ' [min]'
                    else:
                        yvals = (yvals-yvals[0])  # sec
                        ylabelunit = ' [s]'
                elif axesInPlot[1] == 'freq':  # Mhz
                    yvals = yvals/1.e6
                    ylabelunit = ' [MHz]'
            else: 
                yvals = None
                ylabelunit = None

            sf2 = solFetcher(soltab)
            sf2.selection = selection
            # cycle on tables
            titles = []
            dataCube = []
            weightCube = []
            for Ntab, (vals, coord, selection) in enumerate(sf2.getValuesIter(returnAxes=axisInDiff+axisInCol+axesInPlot)):
                dataCube.append([])
                weightCube.append([])

                # set tile
                titles.append('')
                for axis in coord:
                    if axis in axesInFile+axesInPlot+axisInCol: continue
                    titles[Ntab] += axis+':'+str(coord[axis])+' '
                titles[Ntab] = titles[Ntab][:-1] # remove last ' '

                sf3 = solFetcher(soltab)
                sf3.selection = selection
                # cycle on colors
                
                if not ref in sf.getAxisValues('ant') and not ref is None:
                    logging.error('Reference antenna '+ref+' not found. Using: '+sf.getAxisValues('ant')[1])
                    ref = sf.getAxisValues('ant')[1]

                for Ncol, (vals, weight, coord, selection) in enumerate(sf3.getValuesIter(returnAxes=axisInDiff+axesInPlot, weight=True, reference=ref)):
                    dataCube[Ntab].append([])
                    weightCube[Ntab].append([])
        
                    # differential plot
                    if axisInDiff != []:
                        # find ordered list of axis
                        names = [axis for axis in sf.getAxesNames() if axis in axisInDiff+axesInPlot]
                        if axisInDiff[0] not in names:
                            logging.error("Axis to differentiate (%s) not found." % axisInDiff[0])
                            mpm.wait()
                            return 1
                        if len(coord[axisInDiff[0]]) != 2:
                            logging.error("Axis to differentiate (%s) has too many values, only 2 is allowed." % axisInDiff[0])
                            mpm.wait()
                            return 1

                        # find position of interesting axis
                        diff_idx = names.index(axisInDiff[0])
                        # roll to first place
                        vals = np.rollaxis(vals,diff_idx,0)
                        vals = vals[0] - vals[1]
                        weight = np.rollaxis(weight,diff_idx,0)
                        weight = ((weight[0]==1) & (weight[1]==1))
                        del coord[axisInDiff[0]]
 

                    # add tables if required (e.g. phase/tec)
                    for sfAdd in sfsAdd:
                        newCoord = {}
                        for axisName in coord.keys():
                            if axisName in sfAdd.getAxesNames():
                                if type(coord[axisName]) is np.ndarray:
                                    newCoord[axisName] = coord[axisName]
                                else:
                                    newCoord[axisName] = [coord[axisName]] # avoid being interpreted as regexp, faster
                        sfAdd.setSelection(**newCoord)
                        valsAdd = np.squeeze(sfAdd.getValues(retAxesVals=False, weight=False, reference=ref))
                        if sfAdd.getType() == 'clock':
                            valsAdd = 2. * np.pi * valsAdd * newCoord['freq']
                        elif sfAdd.getType() == 'tec':
                            valsAdd = -8.44797245e9 * valsAdd / newCoord['freq']
                        else:
                            logging.warning('Only Clock or TEC can be added to solutions. Ignoring: '+sfAdd.getType()+'.')
                            continue

                        # If clock/tec are single pol then duplicate it (TODO)
                        # There still a problem with commonscalarphase and pol-dependant clock/tec
                        #but there's not easy way to combine them
                        if not 'pol' in sfAdd.getAxesNames() and 'pol' in sf.getAxesNames():
                            # find pol axis positions
                            polAxisPos = sf.getAxesNames().key_idx('pol')
                            # create a new axes for the table to add and duplicate the values
                            valsAdd = np.addaxes(valsAdd, polAxisPos)

                        if valsAdd.shape != vals.shape:
                            logging.error('Cannot combine the table '+sfAdd.getType()+' with '+sf4.getType()+'. Wrong shape.')
                            mpm.wait()
                            return 1

                        vals += valsAdd

                    # normalize
                    if (sf.getType() == 'phase' or sf.getType() == 'scalarphase'):
                        vals = normalize(vals)

                    # unwrap if required
                    if (sf.getType() == 'phase' or sf.getType() == 'scalarphase') and dounwrap:
                        vals = unwrap(vals)
                    
                    # is user requested axis in an order that is different from h5parm, we need to transpose
                    if len(axesInPlot) == 2:
                        if sf3.getAxesNames().index(axesInPlot[0]) < sf3.getAxesNames().index(axesInPlot[1]): vals = vals.T

                    dataCube[Ntab][Ncol] = np.ma.masked_array(vals, mask=(weight == 0))

            # if dataCube too large (> 500 MB) do not go parallel
            if np.array(dataCube).nbytes > 1024*1024*500: 
                logging.debug('Big plot, parallel not possible.')
                plot(Nplots, NColFig, figSize, cmesh, axesInPlot, axisInTable, xvals, yvals, xlabelunit, ylabelunit, datatype, prefix+filename, titles, log, dataCube, minZ, maxZ, plotflag, makeMovie, antCoords, None)
            else:
                mpm.put([Nplots, NColFig, figSize, cmesh, axesInPlot, axisInTable, xvals, yvals, xlabelunit, ylabelunit, datatype, prefix+filename, titles, log, dataCube, minZ, maxZ, plotflag, makeMovie, antCoords])
            if makeMovie: pngs.append(prefix+filename+'.png')

        mpm.wait()

        if makeMovie:
            def long_substr(strings):
                """
                Find longest common substring
                """
                substr = ''
                if len(strings) > 1 and len(strings[0]) > 0:
                    for i in range(len(strings[0])):
                        for j in range(len(strings[0])-i+1):
                            if j > len(substr) and all(strings[0][i:i+j] in x for x in strings):
                                substr = strings[0][i:i+j]
                return substr
            movieName = long_substr(pngs)
            assert movieName != '' # need a common prefix, use prefix keyword in case
            logging.info('Making movie: '+movieName)
            # make every movie last 20 sec, min one second per slide
            fps = np.ceil(len(pngs)/200.)
            ss="mencoder -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:vbitrate=6160000:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:"+\
                    "mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq -mf type=png:fps="+str(fps)+" -nosound -o "+movieName.replace('__tmp__','')+".mpg mf://"+movieName+"*  > mencoder.log 2>&1"
            os.system(ss)
            #print ss
            for png in pngs: os.system('rm '+png)

    return 0
