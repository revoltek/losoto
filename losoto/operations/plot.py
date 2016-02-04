#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implements basic plotting
# WEIGHT: flag-only compliant

import logging
from losoto.operations_lib import *

import multiprocessing
inQueue = multiprocessing.JoinableQueue()

logging.debug('Loading PLOT module.')

class multiThread(multiprocessing.Process):
    """
    This class is a working thread which load parameters from a queue and
    run the flagging on a chunk of data
    """

    def __init__(self, inQueue):
        multiprocessing.Process.__init__(self)
        self.inQueue = inQueue

    def run(self):

        while True:
            parms = self.inQueue.get()

            # poison pill
            if parms is None:
                self.inQueue.task_done()
                break

            self.plot(*parms)
            self.inQueue.task_done()
    
    def plot(self, Nplots, cmesh, axesInPlot, axisInTable, xvals, yvals, xlabelunit, ylabelunit, datatype, filename, titles, log, dataCube, weightCube, minZ, maxZ, plotflag, makeMovie):
        from itertools import cycle, chain
        import numpy as np
        # avoids error if re-setting "agg" a second run of plot
        if not 'matplotlib' in sys.modules:
            import matplotlib as mpl
            mpl.rc('font',size =8 )
            mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
            mpl.use("Agg")
        import matplotlib.pyplot as plt # after setting "Agg" to speed up

        Nr = int(np.ceil(np.sqrt(Nplots)))
        Nc = int(np.ceil(np.float(Nplots)/Nr))
        if makeMovie: figgrid, axa = plt.subplots(Nc, Nr, figsize=(5+2*Nc,4+1*Nr), sharex=True, sharey=True)
        else: figgrid, axa = plt.subplots(Nc, Nr, figsize=(10+3*Nc,8+2*Nr), sharex=True, sharey=True)
        if Nplots == 1: axa = np.array([axa])
        figgrid.subplots_adjust(hspace=0, wspace=0)
        axaiter = chain.from_iterable(axa)

        # axes label 
        if len(axa.shape) == 1: # only one row 
            [ax.set_xlabel(axesInPlot[0]+xlabelunit) for ax in axa[:]]
            if cmesh:
                axa[0].set_ylabel(axesInPlot[1]+ylabelunit)
            else:
                axa[0].set_ylabel(datatype)
        else:
            [ax.set_xlabel(axesInPlot[0]+xlabelunit) for ax in axa[-1,:]]
            if cmesh:
                [ax.set_ylabel(axesInPlot[1]+ylabelunit) for ax in axa[:,0]]
            else:
                [ax.set_ylabel(datatype) for ax in axa[:,0]]

        for Ntab, ax in enumerate(axa.flatten()):
            
            ax.text(.5, .9, titles[Ntab], horizontalalignment='center',fontsize=8,transform=ax.transAxes)
           
            # set log scales if activated
            if 'X' in log: ax.set_xscale('log')
            if 'Y' in log: ax.set_yscale('log')

            # set colors (red reserved for flags)
            colors = cycle(['g', 'b', 'c', 'm', 'y', 'k'])

            for Ncol, data in enumerate(dataCube[Ntab]):

                # set color
                color = next(colors)
                vals = dataCube[Ntab][Ncol]
                weight = weightCube[Ntab][Ncol]

                # plotting
                if cmesh:
                    if minZ == 0: minZ = None
                    if maxZ == 0: maxZ = None
                    if plotflag:
                        vals = np.ma.masked_array(vals, mask=(weight == 0))
                    # if user gives axes names in "wrong" order adapat the values
                    # pcolorfast do not check if x,y,val axes lenghts are coherent
#                    if sf4.getAxesNames().index(axesInPlot[0]) < sf4.getAxesNames().index(axesInPlot[1]): vals = vals.T
                    if log: ax.pcolormesh(xvals, yvals , np.log10(vals), vmin=minZ, vmax=maxZ)
                    else: ax.pcolormesh(xvals, yvals, vals, vmin=minZ, vmax=maxZ)
                    ax.axis([xvals.min(), xvals.max(), yvals.min(), yvals.max()])

                    #plt.colorbar(label=sf.getType())
                else:
                    ax.plot(xvals[np.where(weight!=0)], vals[np.where(weight!=0)], 'o', color=color, markersize=3)
                    if plotflag: ax.plot(xvals[np.where(weight==0)], vals[np.where(weight==0)], 'rx', markersize=3) # plot flagged points
                    if minZ != 0:
                        plt.ylim(ymin=minZ)
                    if maxZ != 0:
                        plt.ylim(ymax=maxZ)

        logging.info("Saving "+filename+'.png')
        if axisInTable != []: plt.savefig(filename+'.png', bbox_inches='tight')
        else: plt.savefig(filename+'.png')
        plt.close()


def run( step, parset, H ):

    import os
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
    minZ, maxZ = parset.getDoubleVector('.'.join(["LoSoTo.Steps", step, "MinMax"]), [0,0] )
    # the axis to plot on one page - e.g. ant to get all antenna's on one plot #
    axisInTable = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "TableAxis"]), [] )
    # the axis to plot in different colours - e.g. pol to get all correlations on one plot #
    axisInCol = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "ColorAxis"]), [] )
    # log='XYZ' to set which axes to put in Log
    log = parset.getString('.'.join(["LoSoTo.Steps", step, "Log"]), "" )
    plotflag = parset.getBool('.'.join(["LoSoTo.Steps", step, "PlotFlag"]), False )
    dounwrap = parset.getBool('.'.join(["LoSoTo.Steps", step, "Unwrap"]), False )
    ref = parset.getString('.'.join(["LoSoTo.Steps", step, "Reference"]), '' )
    tablesToAdd = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Add"]), [] )
    makeMovie = parset.getBool('.'.join(["LoSoTo.Steps", step, "MakeMovie"]), False )
    prefix = parset.getString('.'.join(["LoSoTo.Steps", step, "Prefix"]), '' )
    ncpu = parset.getInt('.'.join(["LoSoTo.Ncpu"]), 1 )

    if makeMovie: 
        prefix = prefix+'__tmp__'

    if os.path.exists(os.path.dirname(prefix)) != '' and not os.path.exists(os.path.dirname(prefix)):
        logging.debug('Creating '+os.path.dirname(prefix)+'.')
        os.makedirs(os.path.dirname(prefix))

    if ref == '': ref = None
    sfsAdd = [ solFetcher(soltab) for soltab in openSoltabs(H, tablesToAdd) ]

    # start processes for multi-thread
    logging.debug('Spowning %i threads...' % ncpu)
    for i in xrange(ncpu):
        t = multiThread(inQueue)
        t.start()

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

        if len(set(axisInTable+axesInPlot+axisInCol)) != len(axisInTable+axesInPlot+axisInCol):
            logging.error('Axis defined multiple times.')
            return 1

        if len(axisInTable) > 1 or len(axisInCol) > 1:
            logging.error('Too many AxisInTable/AxisInCol, they must be at most one each.')
            return 1

        # all axes that are not iterated by anything else
        axesInFile = sf.getAxesNames()
        for axis in axisInTable+axesInPlot+axisInCol:
            axesInFile.remove(axis)
 
        # set subplots scheme
        if axisInTable != []:
            Nplots = sf.getAxisLen(axisInTable[0])
        else:
            Nplots = 1
            
        datatype = sf.getType()

        # cycle on files
        if makeMovie: pngs = [] # store png filenames
        for vals, coord, selection in sf.getValuesIter(returnAxes=axisInTable+axisInCol+axesInPlot):
            
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
                # axis vals (they are always the same, regular arrays)
                yvals = coord[axesInPlot[1]]
                # same as above but for y-axis
                if axesInPlot[1] == 'ant':
                    yvals = np.arange(len(yvals))

                if len(xvals) <= 1 or len(yvals) <=1:
                    logging.error('3D plot must have more then one value per axes.')
                    return 1

                ylabelunit=''
                if axesInPlot[1] == 'time':
                    if len(yvals) > 3600:
                        yvals = (yvals-yvals[0])/3600.  # hrs
                        ylabelunit = ' [hr]'
                    elif len(yvals) > 60:
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
            for Ntab, (vals, coord, selection) in enumerate(sf2.getValuesIter(returnAxes=axisInCol+axesInPlot)):
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

                for Ncol, (vals, weight, coord, selection) in enumerate(sf3.getValuesIter(returnAxes=axesInPlot, weight=True, reference=ref)):
                    dataCube[Ntab].append([])
                    weightCube[Ntab].append([])

                    # add tables if required (e.g. phase/tec)
                    for sfAdd in sfsAdd:
                        newCoord = {}
                        for axisName in coord.keys():
                            if axisName in sfs.getAxesNames():
                                if coord[axisName] is list:
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
                            return 1

                        vals += valsAdd

                    # normalize
                    if (sf.getType() == 'phase' or sf.getType() == 'scalarphase'):
                        vals = normalize(vals)

                    # unwrap if required
                    if (sf.getType() == 'phase' or sf.getType() == 'scalarphase') and dounwrap:
                        vals = unwrap(vals)

                    dataCube[Ntab][Ncol] = vals
                    weightCube[Ntab][Ncol] = weight
    
            inQueue.put([Nplots, cmesh, axesInPlot, axisInTable, xvals, yvals, xlabelunit, ylabelunit, datatype, prefix+filename, titles, log, dataCube, weightCube, minZ, maxZ, plotflag, makeMovie])
            if makeMovie: pngs.append(prefix+filename+'.png')

        # poison pill
        for i in xrange(ncpu):
            inQueue.put(None)
        # finish queue
        inQueue.join()

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
            movieName = long_substr(pngs).replace('__tmp__','')
            assert movieName != '' # need a common prefix, use prefix keyword in case
            logging.info('Making movie: '+movieName)
            # make every movie last 10 sec, min one second per slide
            fps = np.ceil(len(pngs)/20.)
            ss="mencoder -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:vbitrate=6160000:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:"+\
                    "mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq -mf type=png:fps="+str(fps)+" -nosound -o "+movieName+".mpg mf://"+movieName+"*  > mencoder.log 2>&1"
            print ss
            #os.system(ss)
            #for png in pngs: os.system('rm '+png)

    return 0
