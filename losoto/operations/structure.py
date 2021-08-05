#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading STRUCTURE module.')

def _run_parser(soltab, parser, step):
    doUnwrap = parser.getbool( step, 'doUnwrap', False )
    refAnt = parser.getstr( step, 'refAnt', '')
    plotName = parser.getstr( step, 'plotName', '' )
    ndiv = parser.getint( step, 'ndiv', 1 )

    parser.checkSpelling( step, soltab, ['doUnwrap', 'refAnt', 'plotName', 'ndiv'])
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
    if refAnt != '' and refAnt != 'closest' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.warning('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[1]
    if refAnt == '' and doUnwrap:
        logging.error('Unwrap requires reference antenna. Using: '+ants[1])
        refAnt = ants[1]
    if refAnt == '': refAnt = None

    soltab.setSelection(ant='CS*', update=True)

    posAll = soltab.getSolset().getAnt()

    if 'pol' in soltab.getAxesNames(): returnAxes = ['freq','pol','ant','time']
    else: returnAxes = ['freq','ant','time']

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=returnAxes, weight=True, refAnt=refAnt):

        # order positions
        pos = np.array([list(posAll[ant]) for ant in coord['ant']])

        # avg pols
        if 'pol' in soltab.getAxesNames():
            # reorder axes
            vals = reorderAxes( vals, soltab.getAxesNames(), ['pol','ant','freq','time'] )
            weights = reorderAxes( weights, soltab.getAxesNames(), ['pol','ant','freq','time'] )

            vals = np.cos(vals) + 1.j * np.sin(vals)
            vals = np.nansum(vals, axis=0)
            vals = np.angle(vals)
            flags = np.array((weights[0]==0)|(weights[1]==0), dtype=bool)
        else:
            # reorder axes
            vals = reorderAxes( vals, soltab.getAxesNames(), ['ant','freq','time'] )
            weights = reorderAxes( weights, soltab.getAxesNames(), ['ant','freq','time'] )
            flags = np.array((weights==0), dtype=bool)

        # unwrap
        if doUnwrap:
            # remove mean to facilitate unwrapping
            for a, ant in enumerate(coord['ant']):
                if not (flags[a,:,:] == True).all() and ant != refAnt:
                    logging.debug('Unwrapping: '+ant)
                    mean = np.angle( np.nanmean( np.exp(1j*vals[a].flatten()) ))
                    vals[a] -= mean
                    vals[a] = np.mod(vals[a]+np.pi, 2*np.pi) - np.pi
                    vals[a,:,:] = unwrap_2d(vals[a,:,:], flags[a,:,:], coord['freq'], coord['time'])
        
        logging.debug('Computing differential values...')
        t1 = np.ma.array( vals, mask=flags ) # mask flagged data
        dph = t1[np.newaxis]-t1[:,np.newaxis] # ant x ant x freq x time
        D = pos[np.newaxis]-pos[:,np.newaxis] # ant x ant x 3
        D2 = np.triu(np.sqrt(np.sum(D**2,axis=-1))) # calc distance and keep only uppoer triangle larger than 0
        myselect = (D2>0)

        if not doUnwrap:
            logging.debug('Re-normalising...')
            dph = np.mod(dph+np.pi, 2*np.pi) - np.pi
            avgdph = np.ma.average(dph, axis=2) # avg in freq (can do because is between -pi and pi)
            #one extra step to remove most(all) phase wraps, phase wraps disturbe the averaging...
            dph = np.remainder(dph - np.ma.average(avgdph,axis=-1)[:,:,np.newaxis,np.newaxis]+np.pi,2*np.pi) + np.ma.average(avgdph,axis=-1)[:,:,np.newaxis,np.newaxis]-np.pi #center around the avg value
        
        logging.debug('Computing sructure function...')
        avgdph = np.ma.average(dph,axis=2) # avg in freq to reduce noise

        variances = []; pars = []
        avgdph = avgdph[...,avgdph.shape[-1]%ndiv:] # remove a few timeslots to make the array divisible by np.split
        avgdph[avgdph.mask] = np.nan
        for i, avgdphSplit in enumerate( np.split(avgdph, ndiv, axis=-1) ):
            variance = np.nanvar(avgdphSplit, axis=-1)*(np.average(coord['freq'])/150.e6)**2 # get time variance and rescale to 150 MHz

            # linear regression
            #A = np.ones((2,D2[myselect].shape[0]),dtype=float)
            #A[1,:] = np.log10(D2[myselect][~variance.mask])
            #par = np.dot(np.linalg.inv(np.dot(A,A.T)),np.dot(A,np.log10(variance[myselect])))
            mask = np.isnan(variance[myselect])
            A = np.vstack([np.log10(D2[myselect][~mask]), np.ones(len(D2[myselect][~mask]))])
            par = np.linalg.lstsq( A.T, np.log10(variance[myselect][~mask]) )[0] 
            S0 = 10**(-1*par[1]/par[0])
            logging.info(r't%i: beta=%.2f - R_diff=%.3f km' % (i, par[0], S0/1.e3))
            variances.append(variance)
            pars.append(par)

        if plotName != '':
            if plotName.split('.')[-1] != 'png': plotName += '.png' # add png

            if not 'matplotlib' in sys.modules:
                import matplotlib as mpl
                mpl.use("Agg")
            import matplotlib.pyplot as plt
    
            fig = plt.figure()
            fig.subplots_adjust(wspace=0)
            ax = fig.add_subplot(111)
            ax1 = ax.twinx()
            
            for i, variance in enumerate(variances):
                if len(variances) > 1:
                    color = plt.cm.jet(i/float(len(variances)-1)) # from 0 to 1
                else: color = 'black'
                ax.plot(D2[myselect]/1.e3,variance[myselect],marker='o',linestyle='', color=color, markeredgecolor='none', label='T')

                # regression
                par = pars[i]
                x = D2[myselect]
                S0 = 10**(-1*par[1]/par[0])
                if color == 'black': color = 'red' # in case of single color, use red line that is more visible
                ax1.plot(x.flatten()/1.e3, par[0]*np.log10(x.flatten()) + par[1], linestyle='-', color=color, label=r'$\beta=%.2f$ - $R_{\rm diff}=%.2f$ km' % (par[0], S0/1.e3))

            ax.set_xlabel('Distance (km)')
            ax.set_ylabel(r'Phase variance @150 MHz (rad$^2$)')
            ax.set_xscale('log')
            ax.set_yscale('log')

            ymin = np.min(variance[myselect])
            ymax = np.max(variance[myselect])
            ax.set_xlim(xmin=0.1,xmax=3)
            ax.set_ylim(ymin,ymax)
            ax1.set_ylim(np.log10(ymin),np.log10(ymax))
            ax1.legend(loc='lower right', frameon=False)
            ax1.set_yticks([])
        
            logging.warning('Save pic: %s' % plotName)
            plt.savefig(plotName, bbox_inches='tight')

    return 0
