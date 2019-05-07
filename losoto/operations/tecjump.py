#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a jumps remover for TEC solutions
# WEIGHH: flag ready

import logging
from losoto.lib_operations import *

logging.debug('Loading TECJUMP module.')

def _run_parser(soltab, parser, step):
    refAnt = parser.getstr( step, 'refAnt', '' )
    soltabError = parser.getstr( step, 'errorTab', '' )

    parser.checkSpelling( step, soltab, ['refAnt', 'soltabError'])
    return run(soltab, refAnt, soltabError)

def run( soltab, refAnt='', soltabError='' ):
    """
    Remove jumps from TEC solutions.
    WEIGHT: uses the errors.

    Parameters
    ----------
    soltabError : str, optional
        The table name with solution errors. By default it has the same name of soltab with "error" in place of "tec".

    refAnt : str, optional
        Reference antenna for phases. By default None.

    """

    import scipy.ndimage.filters
    import numpy as np
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

    def getPhaseWrapBase(freqs):
        """
        freqs: frequency grid of the data
        return the step size from a local minima (2pi phase wrap) to the others [0]: TEC, [1]: clock
        """
        freqs = np.array(freqs)
        nF = freqs.shape[0]
        A = np.zeros((nF, 2), dtype=np.float)
        A[:, 1] = freqs * 2 * np.pi * 1e-9
        A[:, 0] = -8.44797245e9 / freqs
        steps = np.dot(np.dot(np.linalg.inv(np.dot(A.T, A)), A.T), 2 * np.pi * np.ones((nF, ), dtype=np.float))
        return steps

    class Block(object):
        """
        Implement blocks of contiguous series of good datapoints
        """
        def __init__(self, jump_idx_init, jump_idx_end, vals, vals_e, tec_jump):
            self.tec_jump = tec_jump
            self.jump_idx_init = jump_idx_init
            self.jump_idx_end = jump_idx_end
            self.idx = range(jump_idx_init, jump_idx_end)
            self.idx_exp = range(jump_idx_init-10, jump_idx_end+10)
            self.vals = vals
            self.vals_e = vals_e
            self.len = len(vals)
            self.expand()

        def get_vals(self, idx):
            """
            Return values and errors.
            Idx is already the local index of the block.
            """
            return self.vals_exp[idx], self.vals_e_exp[idx]

        def expand(self):
            """
            predict values+err outside the edges
            TODO: for now it's just nearest
            """
            self.vals_exp = np.concatenate( \
                    ( np.full((10), self.vals[0]) ,
                      self.vals,\
                      np.full((10), self.vals[-1])
                    ) )
            # errors are between 0-1
            self.vals_e_exp = np.concatenate( \
                    ( np.linspace(0.1,1,10) ,
                      self.vals_e,\
                      np.linspace(0.1,1,10)
                    ) )

        def jump(self, jump=0):
            """
            Return values+errors of this block after applying a jump
            """
            self.vals_exp += self.tec_jump*jump

    
    def distance(block1, block2):
        """
        Estimate a distance between two blocks taking advantage of blocks predicted edges
        """
        # find indexes of vals_exp that are common to both blocks
        common_idx_block1 = [i for i, x in enumerate(block1.idx_exp) if x in block2.idx_exp]
        common_idx_block2 = [i for i, x in enumerate(block2.idx_exp) if x in block1.idx_exp]
        if len(common_idx_block1) == 0: return 0 # no overlapping

        v1, v1_e = block1.get_vals(idx=common_idx_block1)
        v2, v2_e = block2.get_vals(idx=common_idx_block2)

        num = np.sum( v1_e * v2_e * np.abs(v1 - v2))
        den = np.sum( v1_e * v2_e )

        return num/den

    def global_dispersion(blocks):
        """
        Calculate a "dispersion" of the time serie summing up the distances between blocks
        """
        dispersion = 0
        for block1, block2 in zip(blocks[:-1], blocks[1:]):
            dispersion += distance(block1, block2)

        return dispersion


    if soltab.getType() != 'tec':
        logging.error('TECJUMP works only on tec solutions.')
        return 1

    if soltabError == '': soltabError = soltab.name.replace('tec','error')
    solset = soltab.getSolset()
    soltab_e = solset.getSoltab(soltabError)
    try:
        seltab_e = solset.getSoltab(soltabError)
    except:
        logging.error('Cannot fine error solution table %s.' % soltabError)
        return 1
    vals_e_all = soltab_e.getValues( retAxesVals = False, weight = False )

    logging.info("Removing TEC jumps from soltab: "+soltab.name)

    ants = soltab.getAxisValues('ant')
    if refAnt != '' and refAnt != 'closest' and not refAnt in soltab.getAxisValues('ant', ignoreSelection = True):
        logging.error('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[0]
    if refAnt == '': refAnt = None

    # Get the theoretical tec jump
    tec_jump_theory = abs(getPhaseWrapBase([42.308e6,42.308e6+23828e6])[0])

    # Find the average jump on all antennas by averaging all jumps found to be withing 0.5 and 2 times the tec_jump_theory
    vals = soltab.getValues( retAxesVals=False, reference=refAnt )
    timeAxis = soltab.getAxesNames().index('time')
    vals = np.swapaxes(vals, 0, timeAxis)
    vals = vals[1,...] - vals[:-1,...] 
    vals = vals[(vals > tec_jump_theory*1) & (vals < tec_jump_theory*1.5)]
    tec_jump = np.median(vals)

    logging.info('TEC jump - theoretical: %.5f TECU - estimated: %.5f TECU' % (tec_jump_theory, tec_jump))

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes='time', weight = True, reference=refAnt):

        # skip all flagged
        if (weights == 0).all(): continue
        # skip reference
        if (vals[(weights != 0)] == 0).all(): continue

        logging.info('Working on ant: %s' % (coord['ant']))

        vals_init = np.copy( vals ) # backup for final check
        vals_e = np.squeeze(vals_e_all[selection])

        # TODO: remove very large jumps (>1 or <-1)

        i = 0
        while True:

            # find blocks
            vals_diff = np.diff(vals)
            vals_diff = np.concatenate(([100], list(vals_diff), [100]))
            jumps_idx = np.where(np.abs(vals_diff) > tec_jump/2.)[0]

            # no more jumps
            if len(jumps_idx) == 2: break

            # add edges and make blocks
            #jumps_idx = [0]+list(jumps_idx)+[len(vals)]
            blocks = []
            for jump_idx_init, jump_idx_end in zip(jumps_idx[:-1],jumps_idx[1:]):
                blocks.append( Block(jump_idx_init, jump_idx_end, \
                        vals[jump_idx_init:jump_idx_end],  vals_e[jump_idx_init:jump_idx_end], tec_jump=tec_jump) )
    
            # cycle on blocks and merge 
            dispersions = []
            for block in blocks:
                block.jump(-1)
                dispersions.append( global_dispersion(blocks) )
                block.jump(+1) # return to normality

                block.jump(+1)
                dispersions.append( global_dispersion(blocks) )
                block.jump(-1) # return to normality
            #print "dispersions:", dispersions
                
            # find best jump
            idx = dispersions.index( min(dispersions) )
            if idx%2 == 0: best_jump = -1
            else: best_jump = +1
            best_block = blocks[idx//2]
            logging.debug('(Cycle %i - #jumps: %i) - Best jump (%i) on block: %i (len %i)' % (i, len(jumps_idx), best_jump, idx//2, best_block.len) )
                
            # recreate vals with the updated block value
            vals[best_block.idx] = best_block.vals + tec_jump*best_jump

            plot = False
            if plot:
                import matplotlib as mpl
                mpl.use("Agg")
                import matplotlib.pyplot as plt
                plt.plot(coord['time'], vals_init, 'k.')
                plt.plot(coord['time'], vals, 'r.')
                plt.savefig('test%03i.png' % i)
                plt.clf()
            i+=1

        # check that this is the closest value to the global minimum
        # (this assumes that the majority of the points are correct)
        distances = []; jumps = []
        for jump in range(-5,5):
            distances.append( np.sum( np.abs(vals_init - (vals + tec_jump*jump) ) ) )
            jumps.append(jump)
        idx = distances.index( min(distances) )
        print('Rescaling all values by %i jumps.' % (jumps[idx]) )
        vals += tec_jump*jumps[idx]

        # set back to 0 the values for flagged data
        vals[weights == 0] = 0
        soltab.setValues(vals, selection)
        soltab.addHistory('TECJUMP')

    return 0


