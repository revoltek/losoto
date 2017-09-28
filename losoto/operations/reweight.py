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


def rolling_window_lastaxis(a, window):
    """Directly taken from Erik Rigtorp's post to numpy-discussion.
    <http://www.mail-archive.com/numpy-discussion@scipy.org/msg29450.html>"""
    import numpy as np

    if window < 1:
       raise ValueError, "`window` must be at least 1."
    if window > a.shape[-1]:
       raise ValueError, "`window` is too long."
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def estimate_weights_window(sindx, vals, nshort, nlong, type, outQueue):
	"""
	Set weights using a median-filter method

	Parameters
	----------
	sindx: int
		Index of station
	vals: array
		Array of values
	nshort: odd int
		Size of short time window
	nlong: odd int
		Size of long time window
	typ: str
		Type of values (e.g., 'phase')

	"""
	import numpy as np

	pad_width = [(0, 0)] * len(vals.shape)
	pad_width[-1] = ((nshort-1)/2, (nshort-1)/2)
	if type == 'phase':
		# Change phase to real/imag (N-1)/2
		real = np.cos(vals)
		imag = np.sin(vals)
		pad_real = np.pad(real, pad_width, 'constant', constant_values=(np.nan,))
		pad_imag = np.pad(imag, pad_width, 'constant', constant_values=(np.nan,))

		# Median smooth with short window and subtract to de-trend
		med_real = np.nanmedian(rolling_window_lastaxis(pad_real, nshort), axis=-1)
		med_imag = np.nanmedian(rolling_window_lastaxis(pad_imag, nshort), axis=-1)
		diff_real = real - med_real
		diff_imag = imag - med_imag
		c = diff_real + diff_imag*1j
	else:
		# Median smooth with short window and subtract to de-trend
		pad_vals = np.pad(vals, pad_width, 'constant', constant_values=(np.nan,))
		med = np.nanmedian(rolling_window_lastaxis(pad_vals, nshort), axis=-1)
		c = vals - med

	# Calculate standard deviation in larger window
	if np.any(c == 0.0):
		c[np.where(np.abs(c) == 0.0)] = np.nan
	pad_width[-1] = ((nlong-1)/2, (nlong-1)/2)
	pad_c = np.pad(c, pad_width, 'constant', constant_values=(np.nan,))
	stddev = np.nanstd(rolling_window_lastaxis(pad_c, nlong), axis=-1)
	if np.any(stddev == 0.0):
		stddev[np.where(stddev == 0.0)] = np.min(stddev[np.where(stddev > 0.0)])

	# Fudge factor to compensate for reduction in scatter due to median
	# and subtraction (i.e., noise remains in median version and is
	# subtracted off)
	corr_factor = 3.0
	w = 1.0 / np.square(corr_factor*stddev)

	outQueue.put([sindx, w])


def run( soltab, method='uniform', weightVal=1., nshort=3, nlong=251,
	soltabImport='', flagBad=False, ncpu=0 ):
	"""
	This operation resets the weight vals

	Parameters
	----------
	method : str, optional
		One of 'uniform' (single value) or 'window' (sliding window in time).
	weightVal : float, optional
		Set weights to this values (0=flagged), by default 1.
	nshort : int, optional
		Short window size in timeslots for 'window' method.
	nlong : int, optional
		Long window size in timeslots for 'window' method.
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
		elif method == 'window':
			if ncpu == 0:
				import multiprocessing
				ncpu = min(8, multiprocessing.cpu_count()) # can use a lot of memory, so don't use too many cores
			if nshort % 2 == 0:
				logging.error('nshort must be odd')
				return 1
			if nlong % 2 == 0:
				logging.error('nlong must be odd')
				return 1

			tindx = soltab.axesNames.index('time')
			antindx = soltab.axesNames.index('ant')
			vals = soltab.val[:].swapaxes(antindx, 0)
			if tindx == 0:
				tindx = antindx
			mpm = multiprocManager(ncpu, estimate_weights_window)
			for sindx, sval in enumerate(vals):
				if np.all(sval == 0.0):
					# skip reference station
					continue
				mpm.put([sindx, sval.swapaxes(tindx-1, -1), nshort, nlong, soltab.getType()])
			mpm.wait()
			weights = np.ones(vals.shape)
			for (sindx, w) in mpm.get():
				weights[sindx, :] = w.swapaxes(-1, tindx-1)
			weights = weights.swapaxes(0, antindx)

			soltab.addHistory('REWEIGHTED using sliding window with nshort={0} '
				'and nlong={1} timeslots'.format(nshort, nlong))
		soltab.setValues(weights, weight=True)

	if flagBad:
		weights = soltab.getValues(weight = True, retAxesVals = False)
		vals = soltab.getValues(retAxesVals = False)
		if soltab.getType() == 'amplitude': weights[np.where(vals == 1)] = 0
		else: weights[np.where(vals == 0)] = 0
		soltab.setValues(weights, weight=True)

	return 0
