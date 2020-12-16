#!/usr/bin/env python
# -*- coding: utf-8 -*-

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading REWEIGHT module.')

def _run_parser(soltab, parser, step):
    mode = parser.getstr( step, 'mode', 'uniform' )
    weightVal = parser.getfloat( step, 'weightVal', 1. )
    nmedian = parser.getint( step, 'nmedian', 3 )
    nstddev = parser.getint( step, 'nstddev', 251 )
    soltabImport = parser.getstr( step, 'soltabImport', '' )
    flagBad = parser.getbool( step, 'flagBad', False )
    ncpu = parser.getint( '_global', 'ncpu', 0 )

    parser.checkSpelling( step, soltab, ['mode', 'weightVal', 'nmedian', 'nstddev', 'soltabImport', 'flagBad'])
    return run(soltab, mode, weightVal, nmedian, nstddev, soltabImport, flagBad, ncpu)


def _rolling_window_lastaxis(a, window):
    """Directly taken from Erik Rigtorp's post to numpy-discussion.
    <http://www.mail-archive.com/numpy-discussion@scipy.org/msg29450.html>"""
    import numpy as np

    if window < 1:
        raise ValueError("`window` must be at least 1.")
    if window > a.shape[-1]:
        raise ValueError("`window` is too long.")
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def _nancircstd(samples, axis=None, is_phase=True):
    """
    Compute the circular standard deviation

    Based on scipy.stats.circstd

    Parameters
    ----------
    imag : array_like
        Input array.
    axis : int, optional
        Axis along which standard deviations are computed.  The default is
        to compute the standard deviation of the flattened array.
    is_phase : bool, optional
        If True, samples are assumed to be phases. If False, they are assumed
        to be either real or imaginary values

    Returns
    -------
    circstd : float
        Circular standard deviation.
    """
    import numpy as np

    if is_phase:
        x1 = np.sin(samples)
        x2 = np.cos(samples)
    else:
        x1 = samples
        x2 = np.sqrt(1.0 - x1**2)
    R = np.hypot(np.nanmean(x1, axis=axis), np.nanmean(x2, axis=axis))

    return np.sqrt(-2*np.log(R))


def _estimate_weights_window(sindx, vals, nmedian, nstddev, stype, outQueue):
    """
    Set weights using a median-filter method

    Parameters
    ----------
    sindx: int
        Index of station
    vals: array
        Array of values
    nmedian: odd int
        Size of median time window
    nstddev: odd int
        Size of stddev time window
    stype: str
        Type of values (e.g., 'phase')

    """
    import numpy as np
    from scipy.ndimage import generic_filter

    pad_width = [(0, 0)] * len(vals.shape)
    pad_width[-1] = (int((nmedian-1)/2), int((nmedian-1)/2))
    if stype == 'phase' or stype == 'rotation':
        # Median smooth and subtract to de-trend
        if nmedian > 0:
            # Convert to real/imag
            real = np.cos(vals)
            pad_real = np.pad(real, pad_width, 'constant', constant_values=(np.nan,))
            med_real = np.nanmedian(_rolling_window_lastaxis(pad_real, nmedian), axis=-1)
            real -= med_real
            real[real < -1.0] = -1.0
            real[real > 1.0] = 1.0

            imag = np.sin(vals)
            pad_imag = np.pad(imag, pad_width, 'constant', constant_values=(np.nan,))
            med_imag = np.nanmedian(_rolling_window_lastaxis(pad_imag, nmedian), axis=-1)
            imag -= med_imag
            imag[imag < -1.0] = -1.0
            imag[imag > 1.0] = 1.0

            # Calculate standard deviations
            pad_width[-1] = (int((nstddev-1)/2), int((nstddev-1)/2))
            pad_real = np.pad(real, pad_width, 'constant', constant_values=(np.nan,))
            stddev1 = _nancircstd(_rolling_window_lastaxis(pad_real, nstddev), axis=-1, is_phase=False)
            pad_imag = np.pad(imag, pad_width, 'constant', constant_values=(np.nan,))
            stddev2 = _nancircstd(_rolling_window_lastaxis(pad_imag, nstddev), axis=-1, is_phase=False)
            stddev = stddev1 + stddev2
        else:
            phase = normalize_phase(vals)

            # Calculate standard deviation
            pad_width[-1] = (int((nstddev-1)/2), int((nstddev-1)/2))
            pad_phase = np.pad(phase, pad_width, 'constant', constant_values=(np.nan,))
            stddev = _nancircstd(_rolling_window_lastaxis(pad_phase, nstddev), axis=-1)
    else:
        if stype == 'amplitude':
            # Assume lognormal distribution for amplitudes
            vals = np.log(vals)

        # Median smooth and subtract to de-trend
        if nmedian > 0:
            pad_vals = np.pad(vals, pad_width, 'constant', constant_values=(np.nan,))
            med = np.nanmedian(_rolling_window_lastaxis(pad_vals, nmedian), axis=-1)
            vals -= med

        # Calculate standard deviation in larger window
        pad_width[-1] = (int((nstddev-1)/2), int((nstddev-1)/2))
        pad_vals = np.pad(vals, pad_width, 'constant', constant_values=(np.nan,))
        stddev = np.nanstd(_rolling_window_lastaxis(pad_vals, nstddev), axis=-1)

    # Check for periods where standard deviation is zero or NaN and replace
    # with min value to prevent inf in the weights. Also limit weights to
    # float16
    zero_scatter_ind = np.where(np.logical_or(np.isnan(stddev), stddev == 0.0))
    if len(zero_scatter_ind[0]) > 0:
        good_ind = np.where(~np.logical_or(np.isnan(stddev), stddev == 0.0))
        stddev[zero_scatter_ind] = np.min(stddev[good_ind])
    if nmedian > 0:
        fudge_factor = 2.0  # factor to compensate for smoothing
    else:
        fudge_factor = 1.0
    w = 1.0 / np.square(stddev*fudge_factor)

    # Rescale to fit in float16
    float16max = 65504.0
    if np.max(w) > float16max:
        w *= float16max / np.max(w)

    outQueue.put([sindx, w])


def run( soltab, mode='uniform', weightVal=1., nmedian=3, nstddev=251,
    soltabImport='', flagBad=False, ncpu=0 ):
    """
    Change the the weight values.

    Parameters
    ----------
    mode : str, optional
        One of 'uniform' (single value), 'window' (sliding window in time), or 'copy' (copy from another table), by default 'uniform'.
    weightVal : float, optional
        Set weights to this values (0=flagged), by default 1.
    nmedian : odd int, optional
        Median window size in number of timeslots for 'window' mode.
        If nonzero, a median-smoothed version of the input values is
        subtracted to detrend them. If 0, no smoothing or subtraction is
        done, by default 3.
    nstddev : odd int, optional
        Standard deviation window size in number of timeslots for 'window' mode, by default 251.
    soltabImport : str, optional
        Name of a soltab. Copy weights from this soltab (must have same axes shape), by default none.
    flagBad : bool, optional
        Re-apply flags to bad values (1 for amp, 0 for other tables), by default False.
    """

    import numpy as np

    logging.info("Reweighting soltab: "+soltab.name)

    if mode == 'copy':
        if soltabImport == '':
            logging.error('In copy mode a soltabImport must be specified.')
            return 1
        solset = soltab.getSolset()
        soltabI = solset.getSoltab(soltabImport)
        soltabI.selection = soltab.selection

        weights, axes = soltab.getValues(weight = True)
        weightsI, axesI = soltabI.getValues(weight = True)
        if list(axes.keys()) != list(axesI.keys()) or weights.shape != weightsI.shape:
            logging.error('Impossible to merge: two tables have with different axes values.')
            return 1
        weightsI[ np.where(weights == 0) ] = 0.
        soltab.setValues(weightsI, weight=True)
        soltab.addHistory('WEIGHT imported from '+soltabI.name+'.')

    elif mode == 'uniform':
        soltab.addHistory('REWEIGHTED to '+str(weightVal)+'.')
        soltab.setValues(weightVal, weight=True)

    elif mode == 'window':
        if nmedian !=0 and nmedian % 2 == 0:
            logging.error('nmedian must be odd')
            return 1
        if nstddev % 2 == 0:
            logging.error('nstddev must be odd')
            return 1

        tindx = soltab.axesNames.index('time')
        antindx = soltab.axesNames.index('ant')
        vals = soltab.val[:].swapaxes(antindx, 0)
        if tindx == 0:
            tindx = antindx
        mpm = multiprocManager(ncpu, _estimate_weights_window)
        for sindx, sval in enumerate(vals):
            if np.all(sval == 0.0) or np.all(np.isnan(sval)):
                # skip reference station
                continue
            mpm.put([sindx, sval.swapaxes(tindx-1, -1), nmedian, nstddev, soltab.getType()])
        mpm.wait()
        weights = np.ones(vals.shape)
        for (sindx, w) in mpm.get():
            weights[sindx, :] = w.swapaxes(-1, tindx-1)
        weights = weights.swapaxes(0, antindx)

        soltab.addHistory('REWEIGHTED using sliding window with nmedian={0} '
            'and nstddev={1} timeslots'.format(nmedian, nstddev))
        soltab.setValues(weights, weight=True)

    if flagBad:
        weights = soltab.getValues(weight = True, retAxesVals = False)
        vals = soltab.getValues(retAxesVals = False)
        if soltab.getType() == 'amplitude': weights[np.where(vals == 1)] = 0
        else: weights[np.where(vals == 0)] = 0
        soltab.setValues(weights, weight=True)

    return 0
