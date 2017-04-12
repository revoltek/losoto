#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a jumps remover for TEC solutions
# WEIGHH: flag ready

import logging
from losoto.operations_lib import *

logging.debug('Loading TECJUMP module.')

def run( step, parset, H ):

    import scipy.ndimage.filters
    import numpy as np
    from losoto.h5parm import solFetcher, solWriter
    from scipy.optimize import minimize
    import itertools
    from scipy.interpolate import griddata
    import scipy.cluster.vq as vq

    def robust_std(data, sigma=3):
        """
        Calculate standard deviation excluding outliers
        ok with masked arrays
        """
        return np.std(data[np.where(np.abs(data) < sigma * np.std(data))])

    def mask_interp(vals, mask, method='nearest'):
        """
        return interpolated values for masked elements
        """
        this_vals = vals.copy()
        #this_vals[mask] = np.interp(np.where(mask)[0], np.where(~mask)[0], vals[~mask])
        #this_vals[mask] = griddata(np.where(~mask)[0], vals[~mask], np.where(mask)[0], method)

        # griddata has nan bug with nearest, I need to use vq
        code, dist = vq.vq(np.where(mask)[0], np.where(~mask)[0])
        this_vals[ np.where(mask)[0] ] = this_vals[code]
        
        return this_vals

    tec_jump_val = 0.019628
    maxsize = 300
    clip = 10 # TECs over these amount of jumps are flagged

    soltabs = getParSoltabs( step, parset, H )

    for soltab in openSoltabs( H, soltabs ):

        logging.info("Removing TEC jumps from soltab: "+soltab._v_name)

        sf = solFetcher(soltab)
        sw = solWriter(soltab) # remember to flush!

        # TODO: check if it's a Tec table

        # axis selection
        userSel = {}
        for axis in sf.getAxesNames():
            userSel[axis] = getParAxis( step, parset, H, axis )
        sf.setSelection(**userSel)

        for vals, weights, coord, selection in sf.getValuesIter(returnAxes='time', weight=True):

            # skip all flagged
            if (weights == 0).all(): continue
            # skip reference
            if (np.diff(vals[(weights == 1)]) == 0).all(): continue

            # kill large values
#            weights[abs(vals/tec_jump_val)>clip] = 0

            # interpolate flagged values to get resonable distances
#            vals = mask_interp(vals, mask=(weights == 0))/tec_jump_val
            # add edges to allow intervals to the borders
#            vals = np.insert(vals, 0, vals[0])
#            vals = np.insert(vals, len(vals), vals[-1])

            vals = np.fmod(vals,tec_jump_val)

#            def find_jumps(d_vals):
#                # jump poistion finder
#                d_smooth = scipy.ndimage.filters.median_filter( mask_interp(d_vals, mask=(abs(d_vals)>0.8)), 21 )
#                d_vals -= d_smooth
#                jumps = list(np.where(np.abs(d_vals) > 1.)[0])
#                return [0]+jumps+[len(d_vals)-1] # add edges
#
#            class Jump(object):
#                def __init__(self, jumps_idx, med):
#                    self.idx_left = jumps_idx[0]
#                    self.idx_right = jumps_idx[1]
#                    self.jump_left = np.rint(d_vals[self.idx_left])
#                    self.jump_right = np.rint(d_vals[self.idx_right])
#                    self.size = self.idx_right-self.idx_left
#                    self.hight = np.median(vals[self.idx_left+1:self.idx_right+1]-med)
#                    if abs((self.hight-self.jump_left)-med) > abs((self.hight-self.jump_right)-med):
#                        self.closejump = self.jump_right
#                    else:
#                        self.closejump = self.jump_left
#    
#            i = 0
#            while i<len(coord['time']):
#                # get tec[i] - tec[i+1], i.e. the derivative assuming constant timesteps
#                # this is in units of tec_jump_val!
#                d_vals = np.diff(vals)
#                # get jumps idx, idx=n means a jump beteen val n and n+1
#                jumps_idx = find_jumps(d_vals)
#
#                # get regions
#                med = np.median(vals)
#                jumps = [Jump(jump_idx, med) for jump_idx in zip( jumps_idx[:-1], jumps_idx[1:] )]
#                jumps = [jump for jump in jumps if jump.closejump != 0]
#                jumps = [jump for jump in jumps if jump.size != 0] # prevent bug on edges
#                jumps = [jump for jump in jumps if jump.size < maxsize]
#
#                jumps.sort(key=lambda x: (np.abs(x.size), x.hight), reverse=False) #smallest first
#                #print [(j.hight, j.closejump) for j in jumps]
#
#                plot = False
#                if plot:
#                    import matplotlib.pyplot as plt
#                    fig, ((ax1, ax2, ax3)) = plt.subplots(3, 1, sharex=True)
#                    fig.subplots_adjust(hspace=0)
#                    d_smooth = scipy.ndimage.filters.median_filter( mask_interp(d_vals, mask=(abs(d_vals)>0.8)), 31 )
#                    ax1.plot(d_vals,'k-')
#                    ax2.plot(d_smooth,'k-')
#                    ax3.plot(vals, 'k-')
#                    [ax3.axvline(jump_idx+0.5, color='r', ls=':') for jump_idx in jumps_idx]
#                    ax1.set_ylabel('d_vals')
#                    ax2.set_ylabel('d_vals - smooth')
#                    ax3.set_ylabel('TEC/jump')
#                    ax3.set_xlabel('timestep')
#                    ax1.set_xlim(xmin=-10, xmax=len(d_smooth)+10)
#                    fig.savefig('plots/%stecjump_debug_%03i' % (coord['ant'], i))
#                i+=1
#
#                if len(jumps) == 0: 
#                    break
#
#                # move down the highest to the side closest to the median
#                j = jumps[0]
#                #print j.idx_left, j.idx_right, j.jump_left, j.jump_right, j.hight, j.closejump
#
#                vals[j.idx_left+1:j.idx_right+1] -= j.closejump
#                logging.debug("%s: Number of jumps left: %i - Removing jump: %i - Size %i" % (coord['ant'], len(jumps_idx)-2, j.closejump, j.size))
                
            # re-create proper vals
#            vals = vals[1:-1]*tec_jump_val
            # set back to 0 the values for flagged data
            vals[weights == 0] = 0

            sw.selection = selection
            sw.setValues(vals)
            sw.setValues(weights, weight=True)

        sw.addHistory('TECJUMP')
        del sf
        del sw
    return 0


