#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a jumps remover for TEC solutions
# WEIGHH: flag ready

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading TECJUMP module.')

def _run_parser(soltab, parser, step):
    refAnt = parser.getstr( step, 'refAnt', '' )
    soltabError = parser.getstr( step, 'errorTab', '' )
    ncpu = parser.getint( '_global', 'ncpu', 0 )

    parser.checkSpelling( step, soltab, ['refAnt', 'errorTab'])
    return run(soltab, refAnt, soltabError, ncpu)


def _run_antenna(vals, vals_e, vals_init, weights, selection, tec_jump, antname, outQueue):
    import itertools, random

    extend = 1 # number of point to eaxtend each block

    class Block(object):
        """
        Implement blocks of contiguous series of good datapoints
        """
        def __init__(self, jump_idx_init, jump_idx_end, vals, vals_e, tec_jump, type='poly'):
            self.id = int(np.mean([jump_idx_init,jump_idx_end])) # not used
            self.tec_jump = tec_jump # empirically calculated expected jump size
            self.jump_idx_init = jump_idx_init # indef of first element of the block
            self.jump_idx_end = jump_idx_end # index of last element of the block
            self.idx = range(jump_idx_init, jump_idx_end) # list of indexes
            self.idx_exp = range(jump_idx_init-extend, jump_idx_end+extend) # add 2 at the edges. These are used to calculate matches with other blocks
            self.vals = vals
            self.vals_e = vals_e
            self.len = len(vals) # lenght of the block

            # fill expected values outside the block boundaries
            self.expand(type)

        def get_vals(self, idx):
            """
            Return values and errors.
            Idx is already the local index of the block.
            """
            return self.vals_exp[idx], self.vals_e_exp[idx]

        def expand(self, type='nearest', order=1, size=4):
            """
            predict values+err outside the edges
            TODO: for now it's just nearest
            """
            if self.len == 1: type='nearest' # if 1 point, then nearest-neighbor
            if (self.len-1) < order: order = self.len-1 # 1st order needts 2 ponts, 2nd order 3 and so on...
            if self.len < size: size = self.len # avoid using more points than available

            if type == 'nearest':
                self.vals_exp = np.concatenate( \
                        ( np.full((extend), self.vals[0]) ,
                          self.vals,\
                          np.full((extend), self.vals[-1])
                        ) )
                # errors are between 0.1 and 1
                #self.vals_e_exp = np.concatenate( \
                #        ( np.linspace(2,0.5,1) ,
                #          self.vals_e,\
                #          np.linspace(0.5,2,1)
                #        ) )
                self.vals_e_exp = np.concatenate( \
                        ( [1.]*extend ,
                          self.vals_e,\
                          [1.]*extend
                        ) )

            elif type == 'poly':
                vals_init = self.vals[:size]
                vals_e_init = self.vals_e[:size]
                idx_init = self.idx[:size]
                vals_end = self.vals[-1*size:]
                vals_e_end = self.vals_e[-1*size:]
                idx_end = self.idx[-1*size:]
                p_init = np.poly1d(np.polyfit(idx_init, vals_init, order, w=1./vals_e_init))
                p_end = np.poly1d(np.polyfit(idx_end, vals_end, order, w=1./vals_e_end))
                self.vals_exp = np.concatenate( \
                        ( p_init(self.idx_exp[:extend]) ,
                          self.vals,\
                          p_end(self.idx_exp[-1*extend:])
                        ) )
                # errors are between 100% and 200%
                self.vals_e_exp = np.concatenate( \
                        ( np.linspace(2,1,extend) ,
                          self.vals_e,\
                          np.linspace(1,2,extend)
                        ) )

        def jump(self, jump=0):
            """
            Return values of this block after applying a jump
            """
            self.vals_exp += self.tec_jump*jump

    
#    def distance(block1, block2):
#        """
#        Estimate a distance between two blocks taking advantage of blocks predicted edges
#        """
#        # find indexes of vals_exp that are common to both blocks
#        common_idx_block1 = [i for i, x in enumerate(block1.idx_exp) if x in block2.idx_exp]
#        common_idx_block2 = [i for i, x in enumerate(block2.idx_exp) if x in block1.idx_exp]
#        if len(common_idx_block1) == 0: return 0 # no overlapping
#
#        v1, v1_e = block1.get_vals(idx=common_idx_block1)
#        v2, v2_e = block2.get_vals(idx=common_idx_block2)
#
#        num = np.sum( v1_e * v2_e * (1+abs(v1 - v2))**(1/3.))
#        den = np.sum( v1_e * v2_e )
#
#        return num/den
#
#    def global_potential_old(blocks):
#        """
#        Calculate a "potential" of the time serie summing up the distances between blocks
#        TODO: rewrite this to go through each point, not each block
#        """
#        potential = 0
#        for block1, block2 in zip(blocks[:-1], blocks[1:]):
#            potential += distance(block1, block2)
#
#        return potential

    def global_potential(blocks, idxs):
        """
        Calculate a "potential" of the time serie going thorugh each point
        """
        potentials = 0
        for idx in idxs:
            vals = []
            vals_w = []
            for block in blocks:
                if idx in block.idx_exp:
                    idx_block = block.idx_exp.index(idx)
                    vals.append(block.vals_exp[idx_block])
                    vals_w.append(1./block.vals_e_exp[idx_block])

            if len(vals) == 1: continue # skip inside block
            potential_n = 0
            potential_d = 0
            for idx1, idx2 in itertools.combinations(range(len(vals)),2):
                potential_n += (vals_w[idx1] * vals_w[idx2] * 1/(vals[idx1] - vals[idx2])**2. )
                potential_d += vals_w[idx1] * vals_w[idx2]

            potentials += potential_n/potential_d
                    
        return potentials

    ###############################################################
    
    for i in range(1000):

        # find blocks
        vals_diff = np.diff(vals)
        vals_diff = np.concatenate(([100], list(vals_diff), [100])) # add edges
        jumps_idx = np.where(np.abs(vals_diff) > tec_jump*(2/3.))[0] # find jumps

        # no more jumps
        if len(jumps_idx) == 2: break

        # make block objects
        blocks = []
        for jump_idx_init, jump_idx_end in zip(jumps_idx[:-1],jumps_idx[1:]):
            blocks.append( Block(jump_idx_init, jump_idx_end, \
                    vals[jump_idx_init:jump_idx_end],  vals_e[jump_idx_init:jump_idx_end], tec_jump=tec_jump) )

        # move the small blocks first
        lens = [block.len for block in blocks]
        if len(set(lens)) > 1:
            max_len = sorted(set(lens))[1]
            best_block = blocks[random.choice( np.where( np.array(lens) <= max_len )[0] )] # random block from the two smallest length
        else:
            best_block = random.choice( blocks ) # random block

        #best_block = random.choice(blocks) # random block
        potentials = []
        jump_range = [-5,-4,-3,-2,-1,1,2,3,4,5]
        for jump_size in jump_range:
            best_block.jump(jump_size)
            idx_touse = range(best_block.idx_exp[0], best_block.idx_exp[extend+1]) + \
                        range(best_block.idx_exp[-1*extend-1], best_block.idx_exp[-1]+1) # check only $extend idx around this block
            #print('block:',best_block.idx_exp,idx_touse)
            #print ('Best block (',best_block.id,') idx:',best_block.idx_exp)
            potentials.append( global_potential(blocks, set(idx_touse) ) )
            best_block.jump(-1*jump_size) # return to normality

        #print "potentials:", potentials
            
        # find best jump
        idx = potentials.index( max(potentials) )
        best_jump = jump_range[idx]
        logging.debug('(Ant %s - Cycle %i - #jumps: %i) - Best jump (%i) on block: %i (len %i)' % (antname, i, len(jumps_idx), best_jump, best_block.id, best_block.len) )
            
        # recreate vals with the updated block value
        vals[best_block.idx] = best_block.vals + tec_jump*best_jump

        plot = False
        if plot:
            print ("Preapare plot")
            best_block.vals_exp += tec_jump*best_jump
            import matplotlib as mpl
            mpl.use("Agg")
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(16, 8))
            ax = fig.add_subplot(111)
            ax.plot(vals_init, 'k.')
            for block in blocks:
                if block is best_block: continue
                #ax.plot(block.idx_exp, block.vals_exp, 'b,')
                ax.errorbar(block.idx_exp[-1:], block.vals_exp[-1:], block.vals_e_exp[-1:]/30., color='blue', ecolor='blue', marker=',', linestyle='')
                ax.errorbar(block.idx_exp[:1], block.vals_exp[:1], block.vals_e_exp[:1]/30., color='blue', ecolor='blue', marker=',', linestyle='')
            #ax.plot(coord['time'], vals, 'r.')
            ax.errorbar(range(len(vals)), vals, vals_e/30., color='green', ecolor='red', marker='.', linestyle='')
            ax.set_xlim(0,len(vals))
            ax.set_ylim(-0.5,0.5)
            fig.savefig('jump_%s_%03i.png' % (antname, i), bbox_inches='tight')
            fig.clf()

    # check that this is the closest value to the global minimum
    # (this assumes that the majority of the points are correct)
    zeros = []; jumps = []
    for jump in range(-100,101):
        #print('TEST jump %i' % jump)
        #print((vals_init - (vals + tec_jump*jump)) )[:10]
        #print( np.where( np.abs(vals_init - (vals + tec_jump*jump)) < 1e-5 )[0] )
        zeros.append( len( np.where( np.abs(vals_init - (vals + tec_jump*jump)) < 1e-5 )[0] ) )
        jumps.append(jump)
    idx = zeros.index( max(zeros) )
    logging.info('%s: Rescaling all values by %i jumps.' % (antname, jumps[idx]) )
    vals += tec_jump*jumps[idx]
    
    outQueue.put([vals, weights, selection])



def run( soltab, refAnt='', soltabError='', ncpu=0 ):
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
    from scipy.interpolate import griddata
    import scipy.cluster.vq as vq

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
        logging.warning('Reference antenna '+refAnt+' not found. Using: '+ants[1])
        refAnt = ants[0]
    if refAnt == '': refAnt = None

    # Get the theoretical tec jump
    tec_jump_theory = abs(getPhaseWrapBase([42.308e6,42.308e6+23828e6])[0])

    # Find the average jump on all antennas by averaging all jumps found to be withing 0.5 and 2 times the tec_jump_theory
    vals = soltab.getValues( retAxesVals=False, refAnt=refAnt )
    timeAxis = soltab.getAxesNames().index('time')
    vals = np.swapaxes(vals, 0, timeAxis)
    vals = vals[1,...] - vals[:-1,...] 
    vals = vals[(vals > tec_jump_theory*1) & (vals < tec_jump_theory*1.5)]
    if len(vals) == 0: 
        logging.info('TEC jump - theoretical: %.5f TECU - NO JUMP FOUND' % (tec_jump_theory))
        return 0
    tec_jump = np.nanmedian(vals)

    logging.info('TEC jump - theoretical: %.5f TECU - estimated: %.5f TECU' % (tec_jump_theory, tec_jump))

    mpm = multiprocManager(ncpu, _run_antenna)

    for vals, weights, coord, selection in soltab.getValuesIter(returnAxes='time', weight = True, refAnt=refAnt):

        # skip all flagged
        if (weights == 0).all(): continue
        # skip reference
        if (vals[(weights != 0)] == 0).all(): continue

        logging.info('Working on ant: %s' % (coord['ant']))

        # 1d linear interp to fill flagged data
        vals[np.where(weights == 0)] = np.interp(np.where(weights == 0)[0], np.where(weights != 0)[0], vals[np.where(weights != 0)])

        vals_init = np.copy( vals ) # backup for final check
        vals_e = np.squeeze(vals_e_all[selection])
        vals_e[np.where(weights == 0)] = 1.

        mpm.put([vals, vals_e, vals_init, weights, selection, tec_jump, coord['ant']])

    mpm.wait()
    for (vals, weights, selection) in mpm.get():
        # set back to 0 the values for flagged data
        vals[weights == 0] = 0
        soltab.setValues(vals, selection)
    soltab.addHistory('TECJUMP')

    return 0
