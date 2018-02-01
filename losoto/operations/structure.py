#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

logging.debug('Loading STRUCTURE module.')

def _run_parser(soltab, parser, step):
    doUnwrap = parser.getbool( step, 'doUnwrap', False )
    refAnt = parser.getstr( step, 'refAnt', '')
    plotName = parser.getstr( step, 'plotName', '' )
    ndiv = parser.getint( step, 'ndiv', 1 )
    return run(soltab, doUnwrap, refAnt, plotName, ndiv)


def run( soltab, doUnwrap=False, refAnt='', plotName='', ndiv=1 ):
    """
    Find the structure function from phase solutions of core stations.

    Parameters
    ----------
    doUnwrap : bool, optional

    refAnt : str, optional
        Reference antenna, by default the first.

    plotName : str, optional
        Plot file name, by default no plot.

    ndiv : int, optional
        

    """
    import numpy as np
    from losoto.lib_unwrap import unwrap, unwrap_2d

    logging.info("Find structure function for soltab: "+soltab.name)

    # input check
    solType = soltab.getType()
    if solType != 'phase':
       logging.warning("Soltab type of "+soltab._v_name+" is of type "+solType+", should be phase.")
       return 1

    ants = soltab.getAxisValues('ant')
    if refAnt != '' and not refAnt in ants:
        logging.error('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[0]
    if refAnt == '': refAnt = ants[0]

    soltab.setSelection(ant='CS*', update=True)

    posAll = soltab.getSolset().getAnt()

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['freq','pol','ant','time'], weight=True, reference=refAnt):

        # reorder axes
        vals = reorderAxes( vals, soltab.getAxesNames(), ['pol','ant','freq','time'] )
        weights = reorderAxes( weights, soltab.getAxesNames(), ['pol','ant','freq','time'] )

        # order positions
        pos = np.array([list(posAll[ant]) for ant in coord['ant']])

        # avg pols
        vals = np.cos(vals) + 1.j * np.sin(vals)
        vals = np.nansum(vals, axis=0)
        vals = np.angle(vals)
        flags = np.array((weights[0]==0)|(weights[1]==0), dtype=bool)

        # unwrap
        if doUnwrap:
            for a, ant in enumerate(coord['ant']):
                if not (flags[a,:,:] == True).all():
                    logging.debug('Unwrapping: '+ant)
                    vals[a,:,:] = unwrap_2d(vals[a,:,:], flags[a,:,:], coord['freq'], coord['time'])
        
        logging.debug('Normilising...')
        t1 = np.ma.array( vals, mask=flags ) # mask flagged data
        dph = t1[np.newaxis]-t1[:,np.newaxis] # ant x ant x freq x time
        D = pos[np.newaxis]-pos[:,np.newaxis] # ant x ant x 3
        D2 = np.triu(np.sqrt(np.sum(D**2,axis=-1))) # calc distance and keep only uppoer triangle larger than 0
        myselect = D2>0
        dph = np.remainder(dph+np.pi,2*np.pi)-np.pi #center around 0
        
        if not doUnwrap:
            logging.debug('Re-normalising...')
            avgdph = np.ma.average(dph, axis=2) # avg in freq (can do because is between -pi and pi)
            #one extra step to remove most(all) phase wraps, phase wraps disturbe the averaging...
            dph = np.remainder(dph - np.ma.average(avgdph,axis=-1)[:,:,np.newaxis,np.newaxis]+np.pi,2*np.pi) + np.ma.average(avgdph,axis=-1)[:,:,np.newaxis,np.newaxis]-np.pi #center around the avg value
        
        logging.debug('Get sructure function...')
        avgdph = np.ma.average(dph,axis=2) # avg in freq to reduce noise

        variances = []
        avgdph = avgdph[...,avgdph.shape[-1]%ndiv:] # remove a few timeslots to make the array divisible by np.split
        for i, avgdphSplit in enumerate( np.split(avgdph, ndiv, axis=-1) ):
            variance = np.ma.var(avgdphSplit, axis=-1)*(np.average(coord['freq'])/150.e6)**2 # get time variance and rescale to 150 MHz

            # linear regression
            A = np.ones((2,D2[myselect].shape[0]),dtype=float)
            A[1,:] = np.log10(D2[myselect])
            par = np.dot(np.linalg.inv(np.dot(A,A.T)),np.dot(A,np.log10(variance[myselect])))
            S0 = 10**(-1*np.array(par)[0]/np.array(par)[1])
            logging.info(r't%i: $\beta=%.2f$ - $R_{diff}=%.2f$ km' % (i, par[1], S0/1.e3))
            variances.append(variance)

        if plotName != '':
            if plotName.split('.')[-1] != 'png': plotName += '.png' # add png

            if not 'matplotlib' in sys.modules:
                import matplotlib as mpl
                mpl.rc('font',size = 20 )
                mpl.rcParams['xtick.labelsize'] = 20
                mpl.use("Agg")
            import matplotlib.pyplot as plt
    
            fig = plt.figure()
            fig.subplots_adjust(wspace=0)
            ax = fig.add_subplot(111)
            
            for i, variance in enumerate(variances):
                if len(variances) > 1:
                    color = plt.cm.jet(i/float(len(variances)-1)) # from 0 to 1
                else: color = 'black'
                ax.plot(D2[myselect]/1.e3,variance[myselect],marker='o',linestyle='',color=color, markeredgecolor='none', label='T')
    
            ax.set_xlabel('Distance (km)')
            ax.set_ylabel(r'Phase variance @150 MHz (rad$^2$)')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(xmin=0.1,xmax=3)
        
            logging.warning('Save pic: %s' % plotName)
            plt.savefig(plotName, bbox_inches='tight')

    return 0
