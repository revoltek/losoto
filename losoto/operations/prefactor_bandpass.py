#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from losoto.lib_operations import *

logging.debug('Loading PREFACTOR_BANDPASS module.')

def _run_parser(soltab, parser, step):
    chanWidth = parser.getstr( step, 'chanWidth')
    BadSBList = parser.getstr( step, 'BadSBList' , '')
    return run(soltab, chanWidth, BadSBList)


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def median_window_filter(ampl, half_window, threshold):
    ampl_tot_copy = np.copy(ampl)
    ndata = len(ampl)
    flags = np.zeros(ndata, dtype=bool)
    sol = np.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]

    #fix oct 2012
    median_array  = scipy.signal.medfilt(sol,half_window*2.-1)

    sol_flag = np.zeros(ndata+2*half_window, dtype=bool)
    sol_flag_val = np.zeros(ndata+2*half_window, dtype=bool)

    for i in range(half_window, half_window + ndata):
        # Compute median of the absolute distance to the median.
        window = sol[i-half_window:i+half_window+1]
        window_flag = sol_flag[i-half_window:i+half_window+1]
        window_masked = window[~window_flag]

        if len(window_masked) < math.sqrt(len(window)):
            # Not enough data to get accurate statistics.
            continue

        median = np.median(window_masked)
        q = 1.4826 * np.median(np.abs(window_masked - median))

        # Flag sample if it is more than 1.4826 * threshold * the
        # median distance away from the median.
        if abs(sol[i] - median) > (threshold * q):
            sol_flag[i] = True

    mask = sol_flag[half_window:half_window + ndata]

    for i in range(len(mask)):
        if mask[i]:
            ampl_tot_copy[i] = median_array[half_window+i] # fixed 2012
    return ampl_tot_copy


def running_median(ampl,half_window) :

    ampl_tot_copy = np.copy(ampl)

    ndata = len(ampl)
    flags = np.zeros(ndata, dtype=bool)
    sol = np.zeros(ndata+2*half_window)
    sol[half_window:half_window+ndata] = ampl
    std = np.zeros(len(ampl))

    for i in range(0, half_window):
        # Mirror at left edge.
        idx = min(ndata-1, half_window-i)
        sol[i] = ampl[idx]

        # Mirror at right edge
        idx = max(0, ndata-2-i)
        sol[ndata+half_window+i] = ampl[idx]

    for i in range(len(ampl)):
        #print i, i+half_window
        std[i] =  np.median(sol[i:i+(2*half_window)])

    return std


def run( soltab, chanWidth, BadSBList = ''):

    import numpy as np
    import scipy
    import scipy.ndimage

    logging.info("Running prefactor_bandpass on: "+soltab.name)
    solset = soltab.getSolset()

    solType = soltab.getType()
    if solType != 'amplitude':
       logging.warning("Soltab type of "+soltab.name+" is: "+solType+" should be amplitude. Ignoring.")
       return 1

    if BadSBList == '':
      bad_sblist   = []
    else:
      bad_sblist = [int(SB) for SB in BadSBList.strip('\"\'').split(';')]

    logging.info("bad SBs: " + str(bad_sblist))

    amplitude_arraytmp = soltab.val[:] # axes are [time, ant, freq, pol]

    logging.info("Shape of amplitudes array: " + str(np.shape(amplitude_arraytmp)))

    nfreqs = len(soltab.freq[:])
    ntimes = len(soltab.time[:])
    nants = len(soltab.ant[:])

    logging.info("Number of antennas: " + str(len(soltab.ant[:])) +  " of frequencies: " + str(nfreqs) + ", and of times: " + str(ntimes))

    subbandHz = 195.3125e3
    if type(chanWidth) is str:
        letters = [1 for s in chanWidth[::-1] if s.isalpha()]
        indx = len(chanWidth) - sum(letters)
        unit = chanWidth[indx:]
        if unit.strip().lower() == 'hz':
            conversion = 1.0
        elif unit.strip().lower() == 'khz':
            conversion = 1e3
        elif unit.strip().lower() == 'mhz':
            conversion = 1e6
        else:
            logging.error("The unit on chanWidth was not understood.")
            raise ValueError("The unit on chanWidth was not understood.")
        chanWidthHz = float(chanWidth[:indx]) * conversion
    else:
        chanWidthHz = chanWidth
    offsetHz = subbandHz / 2.0 - 0.5 * chanWidthHz
    freqmin = np.min(soltab.freq[:]) + offsetHz # central frequency of first subband
    freqmax = np.max(soltab.freq[:]) + offsetHz # central frequency of last subband
    timeidx = np.arange(ntimes)

    SBgrid = np.floor((soltab.freq[:]-np.min(soltab.freq[:]))/subbandHz)

    freqs_new  = np.arange(freqmin, freqmax+100e3, subbandHz)
    amps_array_flagged = np.zeros( (nants, ntimes, len(freqs_new), 2), dtype='float')
    amps_array = np.zeros( (nants, ntimes, len(freqs_new), 2), dtype='float')
    minscale = np.zeros( nants )
    maxscale = np.zeros( nants )

    if len(freqs_new) < 20:
        logging.error("Frequency span is less than 20 subbands! The filtering will not work!")
        logging.error("Please run the calibrator pipeline on the full calibrator bandwidth.")
        raise ValueError("Frequency span is less than 20 subbands! Amplitude filtering will not work!")
        pass

    # remove the bad subbands given by the user
    logging.info("Have " + str(max(SBgrid)) + " subbands.")
    for bad in bad_sblist:
        logging.info('removing subband: ' + str(bad))
        logging.info('removing subband: ' + str(bad))
        pass

    for antenna_id in range(0,len(soltab.ant[:])):
        for time in range(0,len(soltab.time[:])):
            amp_xx_tmp = np.copy(amplitude_arraytmp[time, antenna_id, :, 0])
            amp_yy_tmp = np.copy(amplitude_arraytmp[time, antenna_id, :, 1])
            freq_tmp = soltab.freq[:]
            assert len(amp_xx_tmp[:]) == len(freq_tmp[:])
            mask_xx = np.not_equal(amp_xx_tmp,1.)
            for bad in bad_sblist:
                mask_xx = np.logical_and(SBgrid!=bad,mask_xx)
            if np.sum(mask_xx)>2:
                amps_xx_tointer = amp_xx_tmp[mask_xx]
                freq_xx_tointer = freq_tmp[mask_xx]
                amps_array_flagged[antenna_id,time,:,0] = np.interp(freqs_new,freq_xx_tointer,amps_xx_tointer)
            elif time>0:
                amps_array_flagged[antenna_id,time,:,0] = amps_array_flagged[antenna_id,(time-1),:,0]
            mask_yy = np.not_equal(amp_yy_tmp,1.)
            for bad in bad_sblist:
                mask_yy = np.logical_and(SBgrid!=bad,mask_yy)
            if np.sum(mask_yy)>2:
                amps_yy_tointer = amp_yy_tmp[mask_yy]
                freq_yy_tointer = freq_tmp[mask_yy]
                amps_array_flagged[antenna_id,time,:,1] = np.interp(freqs_new,freq_yy_tointer,amps_yy_tointer)
            elif time>0:
                amps_array_flagged[antenna_id,time,:,1] = amps_array_flagged[antenna_id,(time-1),:,1]

    ampsoutfile = open('calibrator_amplitude_array.txt','w')
    ampsoutfile.write('# Antenna name, Antenna ID, subband, XXamp, YYamp, frequency\n')
    for antenna_id in range(0,len(soltab.ant[:])):
        amp_xx = np.copy(amps_array_flagged[antenna_id,:,:,0])
        amp_yy = np.copy(amps_array_flagged[antenna_id,:,:,1])

        amp_xx = scipy.ndimage.filters.median_filter(amp_xx, (3,3))
        amp_xx = scipy.ndimage.filters.median_filter(amp_xx, (7,1))
        amp_yy = scipy.ndimage.filters.median_filter(amp_yy, (3,3))
        amp_yy = scipy.ndimage.filters.median_filter(amp_yy, (7,1))

        for i in range(0,len(freqs_new)):
            ampsoutfile.write('%s %s %s %s %s %s\n'%(soltab.ant[antenna_id], antenna_id, i, np.median(amp_xx[:,i], axis=0), np.median(amp_yy[:,i], axis=0), freqs_new[i]) )

        for time in range(0,len(soltab.time[:])):
            amps_array[antenna_id,time,:,0] = np.copy(savitzky_golay(amp_xx[time,:], 17, 2))
            amps_array[antenna_id,time,:,1] = np.copy(savitzky_golay(amp_yy[time,:], 17, 2))

        for i in range(0,len(freqs_new)):
            amps_array[antenna_id,:,i,0] = np.median(amps_array[antenna_id,:,i,0])
            amps_array[antenna_id,:,i,1] = np.median(amps_array[antenna_id,:,i,1])
            pass

    try:
        new_soltab = solset.getSoltab('bandpass')
        new_soltab.delete()
    except:
        pass
    new_soltab = solset.makeSoltab(soltype='amplitude', soltabName='bandpass',
                             axesNames=['ant', 'freq', 'pol'], axesVals=[soltab.ant, freqs_new, ['XX','YY']],
                             vals=amps_array[:,0,:,:],
                             weights=np.ones_like(amps_array[:,0,:,:],dtype=np.float16))
    new_soltab.addHistory('CREATE (by PREFACTOR_BANDPASS operation)')

    return 0
