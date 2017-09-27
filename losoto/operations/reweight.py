#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.operations_lib import *
import logging

logging.debug('Loading REWEIGHT module.')

def run_parser(soltab, parser, step):
	weightVal = parser.getfloat( step, 'weightVal', 1. )
	soltabImport = parser.getstr( step, 'soltabImport', '' )
	flagBad = parser.getbool( step, 'flagBad', False )
	return run(soltab, weightVal, soltabImport, flagBad)


def estimate_weights_median_window(soltab, Nstddev):
	"""
	Set weights using a median-filter method
	
	Parameters
	----------
	Nstddev: int
		Size of window in number of timeslots to use to calculate standard deviation

	"""
	import numpy as np

	# Loop over directions
	directions = soltab.dir[:]
	freqs = soltab.freq[:]
	stations = soltab.ant[:]
	for d in directions:
		for f in freqs:
			for s in stations:
				soltab.setSelection(dir=d, freq=f, ant=s)
				phases = np.squeeze(soltab.val)
				if np.all(phases == 0.0):
					# If reference station, skip
					continue
		
				# Change phase to real/imag
				real = np.cos(phases)
				imag = np.sin(phases)

				# Smooth with short window and subtract
				N = 3
				idx = np.arange(N) + np.arange(len(phases)-N+1)[:, None]
				med = np.median(real[idx], axis=1)
				med_real = np.zeros(phases.shape)
				med_real[0:(N-1)/2] = med[0]
				med_real[(N-1)/2:-(N-1)/2] = med
				med_real[-(N-1)/2:] = med[-1]
				med = np.median(imag[idx], axis=1)
				med_imag = np.zeros(phases.shape)
				med_imag[0:(N-1)/2] = med[0]
				med_imag[(N-1)/2:-(N-1)/2] = med
				med_imag[-(N-1)/2:] = med[-1]
				diff_real = real - med_real
				diff_imag = imag - med_imag
				c = diff_real + diff_imag*1j

				# Calculate standard deviation in larger window
				N = Nstddev
				idx = np.arange(N) + np.arange(len(phases)-N+1)[:, None]
				mad = np.std(c[idx], axis=1)
				mad_c = np.zeros(phases.shape)
				mad_c[0:(N-1)/2] = mad[0]
				mad_c[(N-1)/2:-(N-1)/2] = mad
				mad_c[-(N-1)/2:] = mad[-1]
				mad_c[mad_c == 0.0] = np.min(mad_c[mad_c > 0.0])
				nzeros = np.where(c == 0.0)[0].shape[0]
				if nzeros > 0:
					corr_factor = 2.0 * float(phases.shape[0]) / float(nzeros) # compensate for reduction in scatter due to smoothing
				else:
					corr_factor = 1.0
				weights = 1.0 / np.square(corr_factor*mad_c)
				soltab.setValues(weights, weight=True)

	return soltab
	

def run( soltab, method='uniform', weightVal=1., nwindow=250, soltabImport='', flagBad=False ):
	"""
	This operation reset the weight vals

	Parameters
	----------
	method : str, optional
		One of 'uniform' or 'window'.
	weightVal : float, optional
		Set weights to this values (0=flagged), by default 1.
	nwindow : int, optional
		Window size in timeslots for 'window' method.
	soltabImport : str, optional
		Name of a soltab. Copy weights from this soltab, by default do not copy.
	flagBad : bool, optional
		Re-apply flags to bad values, by default False.
	"""

	import numpy as np

	logging.info("Reweighting soltab: "+soltab.name)

	if soltabImport != '':
		solset = soltab.getSolset()
		soltabI = solset.getSoltab(soltabImport)
		soltabI.selection = soltab.selection

		weights, axes = soltab.getValues(weight = True)
		weightsI, axesI = soltabI.getValues(weight = True)
		if axes.keys() != axesI.keys() or weights.shape != weightsI.shape:
			logging.error('Impossible to merge two tables with different axes values')
			return 1
		weightsI[ np.where(weights == 0) ] = 0.
		soltab.setValues(weightsI, weight=True)
		soltab.addHistory('WEIGHT imported from '+soltabI.name+'.')
	else:
		if method == 'uniform':
			weights = weightVal
			soltab.addHistory('REWEIGHTED to '+str(weightVal)+'.')
			soltab.setValues(weights, weight=True)
		elif method == 'window':
			if soltab.getType() != 'phase':
				logging.error('Median-window weighting only works for phases')
				return 1
			soltab = estimate_weights_median_window(soltab, nwindow)
			soltab.addHistory('REWEIGHTED using sliding window of size {} timeslots'.format(nwindow))

	if flagBad:
		weights = soltab.getValues(weight = True, retAxesVals = False)
		vals = soltab.getValues(retAxesVals = False)
		if soltab.getType() == 'amplitude': weights[np.where(vals == 1)] = 0
		else: weights[np.where(vals == 0)] = 0
		soltab.setValues(weights, weight=True)

	return 0
