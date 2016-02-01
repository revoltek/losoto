#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import numpy.ma as ma
import sys
import logging
from lofar.expion import baselinefitting as fitting


# from pylab import *

def ClockTECfunc(xarray, par):
    delay = np.array([par[1] * 1e-9]).flatten()  # in ns, array has dimension 1, even if scalar
    delayfact = 2 * np.pi * delay[:, np.newaxis] * xarray
    TEC = np.array([par[0]]).flatten()  # dTEC in TECU
    drefract = -8.4479745e9 * TEC[:, np.newaxis] / xarray
    if len(par) > 2:
        return drefract[:, np.newaxis, :] + delayfact[np.newaxis] + par[2]  # returns nTEC x nClock x nFreq

    return drefract[:, np.newaxis, :] + delayfact


def ClockTECfuncAllStations(xarray, par):
    delay = np.array([par[1] * 1e-9]).flatten()  # in ns, array has dimension 1, even if scalar
    delayfact = 2 * np.pi * delay[:, np.newaxis] * xarray
    # print "delayfact",delayfact.shape
    TEC = np.array([par[0]]).flatten()  # dTEC in TECU
    drefract = -8.4479745e9 * TEC[:, np.newaxis] / xarray
    # print "refract",drefract.shape
    if len(par) > 2:
        return drefract + delayfact + par[2][:, np.newaxis]  # returns nTEC x nClock x nFreq

    return drefract + delayfact


def getInitClock(data, freq):
    nF = freq.shape[0]
    avgdata = np.ma.sum(np.cos(data) + 1.j * np.sin(data), axis=0).swapaxes(0, -2)
    avgdata = np.ma.arctan2(np.imag(avgdata), np.real(avgdata))
    nSt = avgdata.shape[0]
    npol = avgdata.shape[2]
    for ist in xrange(nSt):
        for pol in xrange(npol):
            mymask = avgdata[ist, :, pol].mask
            if not hasattr(mymask, '__len__'):
                mymask = np.ones(avgdata[ist, :, pol].shape, dtype=bool) * mymask
            if avgdata[ist, :, pol].count() < 1:
                avgdata[ist, :, pol].mask[0] = False
            # logging.info("mask station %d pol %d "%(ist,pol) +str(mymask))
            # logging.info("average data station %d pol %d "%(ist,pol) +str(avgdata[ist,:,pol]))
            avgdata[ist, :, pol][~mymask] = np.float32(np.unwrap(avgdata[ist, :, pol][~mymask]))
            # logging.info("average unwrapped data station %d pol %d "%(ist,pol) +str(avgdata[ist,:,pol]))
            # logging.info("remainder " +str(np.remainder(avgdata[ist,:,pol]+np.pi,2*np.pi)-np.pi))
    A = np.ones((nF, 2), dtype=np.float)
    A[:, 1] = freq * 2 * np.pi * 1e-9
    return np.ma.dot(np.linalg.inv(np.dot(A.T, A)), np.ma.dot(A.T, avgdata).swapaxes(0, -2))


def getInitPar(
    data,
    dTECArray,
    dClockArray,
    freqs,
    ff=ClockTECfunc,
    ):
    '''initialize paramaters and unwraps data for fit'''

    if np.min(np.absolute(dTECArray)) < 1e-5:
        dTECArray += 0.0001  # just to prevent 0 to be there, since otherwise the fit might get stuck in 0
    nT = dTECArray.shape[0]
    nD = dClockArray.shape[0]
    par = [dTECArray, dClockArray, 0]
    # first check all unwrapping possibilities
    bigdata = ff(freqs, par)  # returns array of shape nT,nD,nF
    wraps = np.ma.around(np.divide(bigdata - data[np.newaxis, np.newaxis], 2 * np.pi))
    difference = bigdata - data[np.newaxis, np.newaxis] - wraps * 2 * np.pi
    index = np.unravel_index(np.ma.argmin(np.ma.var(difference, axis=2)), (nT, nD))  # to find abs and remove offset
    OffsetIn = -1 * np.ma.mean(difference[index])
    par = [dTECArray[index[0]], dClockArray[index[1]], OffsetIn]
    estimate = ff(freqs, par).flatten()
    # logging.info("estimate "+str(estimate))
    wraps = np.ma.around(np.divide(estimate - data, 2 * np.pi))
    data[:] = np.add(2 * np.pi * wraps, data)  # update the data
    return par


def getClockTECFit(
    ph,
    freq,
    stations,
    initSol=[],
    returnResiduals=True,
    chi2cut=1e8,
    ):
    stepFraction = 0.1  # step fraction of a 2pi step for brutforce
    # TODO: add a nonlinear fit after the first guess to arrive to the bottom of the local minima
    nT = ph.shape[0]
    nF = freq.shape[0]
    nSt = ph.shape[2]
    data = ph
    # logging.info("fitting masked data "+str(ph.count(axis=0)))
    tecarray = np.zeros((nT, nSt), dtype=np.float32)
    clockarray = np.zeros((nT, nSt), dtype=np.float32)
    if returnResiduals:
        residualarray = np.zeros((nT, nF, nSt), dtype=np.float32)
    A = np.ones((nF, 2), dtype=np.float)
    A[:, 1] = freq * 2 * np.pi * 1e-9
    A[:, 0] = -8.44797245e9 / freq
    (_base, steps) = getPhaseWrapBase(freq)
    stepdTEC = np.abs(steps[0]) * stepFraction
    stepDelay = np.abs(steps[1]) * stepFraction
    succes = False
    initprevsol = np.zeros(nSt, dtype=bool)
    nrFail = np.zeros(nSt, dtype=int)  # number of fails after last success
    sol = np.zeros((nSt, 2), dtype=np.float)
    #
    prevsol = np.copy(sol)
    for itm in xrange(nT):

        # if itm%100==0 and itm>0:
            # sys.stdout.write(str(itm)+'... '+str(sol[-1,0]-sol[0,0])+' '+str(sol[-1,1]-sol[0,1])+' '+str(sol[-1,-1]-sol[0,-1])+' ')
            # sys.stdout.flush()

        if itm == 0 or not succes:
            for ist in xrange(nSt):
                # very first step
                if itm == 0 or not initprevsol[ist]:
                    if hasattr(initSol, '__len__') and len(initSol) > ist:
                        iTEC1 = initSol[ist, 0]
                        iTEC2 = initSol[ist, 0] + stepdTEC
                        iD1 = initSol[ist, 1]
                        iD2 = initSol[ist, 1] + stepDelay
                    else:
                        if 'CS' in stations[ist]:
                            #  TODO: here it assumes that ref station is a core station
                            iTEC1 = -0.2
                            iTEC2 = 0.2
                            iD1 = -20
                            iD2 = 20
                        else:
                            if 'RS' in stations[ist]:
                                iD1 = -250
                                iD2 = 250
                                iTEC1 = -1
                                iTEC2 = 1
                            else:
                                # large TEC variation for EU stations
                                iD1 = -250
                                iD2 = 250
                                iTEC1 = -5
                                iTEC2 = 5
                            if 'LBA' in stations[ist]:
                                # no init clock possible due to large TEC effect
                                iD1 = -400
                                iD2 = 400

                    logging.info('First %f %f %f %f %f %f ' % (
                        iTEC1,
                        iTEC2,
                        stepdTEC,
                        iD1,
                        iD2,
                        stepDelay,
                        ))
                else:
                    # further steps with non success
                    sol[ist, :] = prevsol[ist, :]
                    iTEC1 = prevsol[ist, 0] - min(1.5, stepdTEC * int(nrFail[ist] / 1))  # /stepFraction
                    iTEC2 = prevsol[ist, 0] + min(1.5, stepdTEC * (int(nrFail[ist] / 1) + 1))  # /stepFraction
                    if not 'CS' in stations[ist]:
                        iD1 = prevsol[ist, 1] - min(100, stepDelay * int(nrFail[ist] / 20))  # /stepFraction
                        iD2 = prevsol[ist, 1] + min(100, stepDelay * (int(nrFail[ist] / 20) + 1))  # /stepFraction
                    else:
                        iD1 = prevsol[ist, 1]
                        iD2 = prevsol[ist, 1] + stepDelay

                # logging.info("Failure %d : %f %f %f %f %f %f %d "%(ist,iTEC1,iTEC2,stepdTEC,iD1,iD2,stepDelay,nrFail)+str(prevsol[ist]))
                dTECArray = np.arange(iTEC1, iTEC2, stepdTEC)
                dClockArray = np.arange(iD1, iD2, stepDelay)
                datatmp = ph[itm, :, ist]
                # logging.info("getting init par for station %d"%ist)
                if datatmp.count() / float(nF) > 0.7:
                    # do brutforce and update data
                    par = getInitPar(datatmp, dTECArray, dClockArray, freq, ClockTECfunc)
                    sol[ist, :] = par[:2]

        # every time slot
        for nr_iter in xrange(2):
            # first do a brute unwrap by commparing with estimated phases
            estimate = ClockTECfuncAllStations(freq, sol.T).reshape((nSt, nF)).T
            wraps = np.ma.around(np.divide(estimate - data[itm], 2 * np.pi))
            # find 2pi jumps in unwrapped data (versus frequency)
            wrapjumps = data[itm, 1:] + 2 * np.pi * wraps[1:] - data[itm, :-1] - 2 * np.pi * wraps[:-1]
            jumps = np.int32(np.logical_and(np.absolute(wrapjumps[1:-1]) > 1.5 * np.pi, np.logical_and(np.absolute(wrapjumps[2:]) < 0.5 * np.pi, np.absolute(wrapjumps[:-2]) < 0.5 * np.pi)))
            jumps *= np.sign(wrapjumps[1:-1]).astype(jumps.dtype)
            # jumps should be subtracted cummulative
            jumps = np.cumsum(jumps, axis=0)
            wraps[2:-1] -= jumps
            wraps[-1] -= jumps[-1]

#            if doplot:
#                for ist in range(nSt):
#                    if np.sum(np.absolute(jumps[:,ist]))>0:
#                        subplot(2,1,2)
#                        plot(data[itm,:,ist]+2*np.pi*wraps[:,ist])
#                        title("iter "+str(nr_iter)+" "+str(ist))
#                        break
#                show()

            data[itm, :] = np.ma.add(2 * np.pi * wraps, data[itm])
            # subtract any "global" jump artifically intrduced
            data[itm, :] += np.ma.around(np.ma.average(estimate - data[itm], axis=0) / (2 * np.pi))[np.newaxis] * 2 * np.pi

            wrapflags = np.ma.absolute(estimate - data[itm, :]) < 4. / (2 * (nr_iter + 1)) * np.pi  # flagging bad freq chans of this time slot per station

            for ist in xrange(nSt):

                if data[itm, :, ist][wrapflags[:, ist]].count() / float(nF) < 0.5:
                    # logging.info("too many data points flagged t=%d st=%d flags=%d wrappedflags=%d"%(itm,ist,data[itm,:,ist].count(),data[itm,:,ist][wrapflags[:,ist]].count()) + str(sol[ist])+" "+str(nr_iter)+" "+str(np.ma.absolute(estimate-data[itm,:])[:ist]))
                    sol[ist] = [-10., -10.]
                    continue

                # adding mask to A
                B = np.ma.array(A[wrapflags[:, ist]], mask=np.concatenate((data[itm, :, ist][wrapflags[:, ist]].mask[:, np.newaxis], data[itm, :, ist][wrapflags[:, ist]].mask[:, np.newaxis]), axis=1))
                sol[ist] = np.ma.dot(np.linalg.inv(np.ma.dot(B.T, B)), np.ma.dot(B.T, data[itm, :, ist][wrapflags[:, ist]])).T

                # check if data has jumps compared with previous timeslot
                if initprevsol[ist] and np.abs((sol[ist, 1] - prevsol[ist, 1]) / steps[1]) > 0.5 and ('LBA' in stations[ist] or np.abs((sol[ist, 0] - prevsol[ist, 0]) / steps[0]) > 0.5):
                    # logging.info("removing jumps %f %f %f"%(sol[ist,1],prevsol[ist,1],steps[1])+" "+str('LBA' in stations[ist])+" "+stations[ist]+str(sol[ist])+" "+str(itm))
                    sol[ist, :] -= np.round((sol[ist, 1] - prevsol[ist, 1]) / steps[1]) * steps

        # calculate chi2 per station
        residual = data[itm] - np.dot(A, sol.T)
        tmpresid = residual - residual[:, 0][:, np.newaxis]  # residuals relative to station 0
        residual = np.ma.remainder(tmpresid + np.pi, 2 * np.pi) - np.pi
        chi2 = np.ma.sum(np.square(np.degrees(residual)), axis=0) / nF

        if returnResiduals:
            residualarray[itm] = residual

        chi2select = np.logical_or(np.array(chi2 > chi2cut), sol[:, 0] < -5) # select bad points
        if np.any(chi2select):
            logging.debug('high chi2 of fit, itm: %d %d ' % (itm, np.sum(chi2select)) + str(sol[chi2select]) + 'stations:' + str(np.arange(nSt)[chi2select]) + ' chi2 ' + str(chi2[chi2select]))
            succes = False
            nrFail[chi2select] += 1
            nrFail[~chi2select] = 0
            prevsol[~chi2select][prevsol[~chi2select] == 0] = sol[~chi2select][prevsol[~chi2select] == 0]  # compensate missing prevsol at first rounds
            prevsol[~chi2select] = 0.5 * prevsol[~chi2select] + 0.5 * sol[~chi2select]  # init solution to 0.5 * this solution + 0.5 previous solution
            initprevsol[~chi2select] = True # once is True it never becomes False
        else:
            succes = True
            prevsol[prevsol == 0] = sol[prevsol == 0]  # compensate missing prevsol at first rounds
            prevsol = 0.5 * prevsol + 0.5 * np.copy(sol)
            initprevsol = np.ones(nSt, dtype=bool)
            nrFail = np.zeros(sol.shape[0], dtype=int)

#        nrFail[chi2select] += 1
#        nrFail[~chi2select] = 0
        # compensate missing prevsol at first rounds
#        prevsol[~chi2select][prevsol[~chi2select] == 0] = sol[~chi2select][prevsol[~chi2select] == 0]
#        prevsol[~chi2select] = 0.5 * prevsol[~chi2select] + 0.5 * sol[~chi2select]  # init solution to 0.5 * this solution + 0.5 * previous solution
#        prevsol[~chi2select] = sol[~chi2select] + sol[~chi2select] - prevsol[~chi2select]  # init solution using the extrapolated value (linear regression assuming x2-x1 == 1)
        tecarray[itm] = sol[:, 0]
        clockarray[itm] = sol[:, 1]
    if returnResiduals:
        return (tecarray, clockarray, residualarray)
    return (tecarray, clockarray)


def getPhaseWrapBase(freqs):
    """
    freqs: frequency grid of the data
    return the step size from a local minima (2pi phase wrap) to the others [0]: TEC, [1]: clock
    """

    nF = freqs.shape[0]
    A = np.zeros((nF, 2), dtype=np.float)
    A[:, 1] = freqs * 2 * np.pi * 1e-9
    A[:, 0] = -8.44797245e9 / freqs
    steps = np.dot(np.dot(np.linalg.inv(np.dot(A.T, A)), A.T), 2 * np.pi * np.ones((nF, ), dtype=np.float))
    basef = np.dot(A, steps) - 2 * np.pi
    return (basef, steps)


def getResidualPhaseWraps(avgResiduals, freqs):
    flags = avgResiduals[:, 10] == 0.
    nSt = avgResiduals.shape[1]
    nF = freqs.shape[0]
    wraps = np.zeros((nSt, ), dtype=np.float)
    tmpflags = flags
    tmpfreqs = freqs[np.logical_not(tmpflags)]
    (tmpbasef, steps) = getPhaseWrapBase(tmpfreqs)
    basef = np.zeros(freqs.shape)
    basef[np.logical_not(tmpflags)] = tmpbasef
    basef = basef.reshape((-1, 1))
    data = avgResiduals[:, :]
    wraps = fitting.fit(data, basef, wraps, flags).flatten()
    return (wraps, steps)


def correctWraps(
    tecarray,
    residualarray,
    freq,
    pos,
    ):
    nT = tecarray.shape[0]
    nSt = tecarray.shape[1]
    flags = tecarray < -5
    resflags = np.logical_or(flags[:, np.newaxis], residualarray == 0)
    maskedresiduals = np.ma.array(residualarray, mask=resflags)
    # avgResiduals=np.average(residualarray,axis=0)
    avgResiduals = np.ma.average(maskedresiduals, axis=0)
    (wraps, steps) = getResidualPhaseWraps(avgResiduals, freq)  # fitting of the wraps from the time-avg residuals
    # step[0] is the step in TEC corresponding to 1 phase wrap and step[1] is in ns (clock)
    wraps = np.round(wraps - wraps[0])  # reference to station 0
    logging.info('wraps from residuals: ' + str(wraps))
    lats = np.degrees(np.arctan2(pos[:, 2], np.sqrt(pos[:, 0] * pos[:, 0] + pos[:, 1] * pos[:, 1])))
    lats -= lats[0]
    lons = np.degrees(np.arctan2(pos[:, 1], pos[:, 0]))
    lons -= lons[0]
    lonlat = np.concatenate((lons, lats)).reshape((2, ) + lons.shape)
    # refine (is it needed for LBA? check results)
    for nr_iter in xrange(2):
        # recreate best TEC at the moment
        TEC = tecarray - tecarray[:, [0]] + steps[0] * (np.round(wraps) - np.round(wraps[0]))
        TEC = np.ma.array(TEC, mask=flags)
        # fit 2d linear TEC screen over stations
        slope = np.ma.dot(np.linalg.inv(np.dot(lonlat, lonlat.T)), np.ma.dot(lonlat, TEC.T))
        # flag bad time steps maybe because TEC screen is not linear
        chi2 = np.ma.sum(np.square(TEC - np.ma.dot(lonlat.T, slope).T), axis=1) / nSt
        chi2select = chi2 < np.ma.average(chi2)
        chi2select = chi2 < np.ma.average(chi2[chi2select])
        # calculate offset per station wrt time-averaged TEC screen
        offsets = -1 * np.ma.average(TEC[chi2select] - np.ma.dot(slope.T, lonlat)[chi2select], axis=0) * 2. * np.pi / steps[0]
        remainingwraps = np.round(offsets / (2 * np.pi))  # -np.round(wraps[stationIndices])
        logging.info('offsets: ' + str(offsets))
        logging.info('avgTEC: ' + str(np.ma.average(TEC[chi2select], axis=0)))
        logging.info('remaining: ' + str(remainingwraps))
        wraps += remainingwraps
        # TODO: remove also the offset before the second cycle
        if np.sum(np.absolute(remainingwraps)) == 0:
            break
    return (offsets, wraps, steps)


def doFit(
    phases,
    mask,
    freqs,
    stations,
    station_positions,
    axes,
    refstIdx='superterp',
    stationSelect='BA',
    flagBadChannels=True,
    flagcut=1.5,
    chi2cut=30000.,
    removePhaseWraps=True,
    combine_pol=False,
    initSol=[],
    ):
    # make sure order of axes is as expected
    stidx = axes.index('ant')
    freqidx = axes.index('freq')
    timeidx = axes.index('time')
    polidx = axes.index('pol')
    data = ma.array(phases, mask=mask).transpose((timeidx, freqidx, stidx, polidx))
    nT = data.shape[0]
    nF = data.shape[1]
    nSt = data.shape[2]
    npol = data.shape[3]
    if npol == 4:
        data = data[:, :, :, (0, 3)]
        npol = 2
    if refstIdx == 'superterp':
        superterpstations = [i for i in stations if i[:5] in [
            'CS002',
            'CS003',
            'CS004',
            'CS005',
            'CS006',
            'CS007',
            ]]
        refstIdx = [i for (i, j) in enumerate(stations) if j in superterpstations]
    if not hasattr(refstIdx, '__len__'):
        refstIdx = [refstIdx]

    # get phases from reference stations
    refdata = np.ma.sum(np.cos(data[:, :, refstIdx, :]) + 1.j * np.sin(data[:, :, refstIdx, :]), axis=2)
    refdata = np.ma.arctan2(np.imag(refdata), np.real(refdata))
    # unwrap around mean
    mymean = np.ma.average(refdata, axis=0)
    refdata = np.ma.remainder(refdata - mymean[np.newaxis] + np.pi, 2 * np.pi) + mymean[np.newaxis] - np.pi
    data -= refdata[:, :, np.newaxis]
    # flag bad channels - test if really nedded if data already flagged
    indices = np.arange(nF)
    if flagBadChannels:
        freqselect = np.ones((nF, ), dtype=np.bool)
        for nr_iter in xrange(2):
            rms = np.ma.std(np.ma.std(refdata, axis=0), axis=1)
            freqselect = rms < flagcut * np.average(rms)
            logging.info('iter %d: flagging %d channels' % (nr_iter, np.sum(np.logical_not(freqselect))))
            freqs = freqs[freqselect]
            data = data[:, freqselect]
            refdata = refdata[:, freqselect]
            logging.info('flagging: ' + str(indices[np.logical_not(freqselect)]))
            indices = indices[freqselect]
        nF = data.shape[1]
    # select stations - can be removed
    if isinstance(stationSelect, str):
        selectstations = [st for st in stations if stationSelect in st]
    else:
        selectstations = list(stations[stationSelect])
    logging.info('%d selected stations: ' % len(selectstations) + str(selectstations))
    stationIndices = np.array([idxst in selectstations for idxst in stations])
    data = data[:, :, stationIndices]
    stations = stations[stationIndices]
    station_positions = station_positions[stationIndices]
    RSstations = [i for (i, j) in enumerate(stations) if 'RS' in j]
    CSstations = [i for (i, j) in enumerate(stations) if 'CS' in j]
    otherstations = [i for (i, j) in enumerate(stations) if not 'RS' in j and not 'CS' in j]
    logging.info('station indices: ' + str(stationIndices) + ' RS ' + str(RSstations))
    nSt = data.shape[2]
    logging.info('%d selected stations: ' % nSt + str(stations))
    # combine polarizationsif requested - needed in HBA
    if combine_pol:
        if npol == 2:
            cdata = np.ma.cos(data) + 1.j * np.ma.sin(data)
            data = np.ma.sum(cdata, axis=3).reshape((nT, nF, nSt, 1))
            data = np.ma.arctan2(np.imag(data), np.real(data))  # np.angle doesnot yet return masked array!!
            npol = 1
    # guess clock, remove from data
    # not in LBA because TEC dominant
    if not 'LBA' in stations[0] and len(initSol) < 1:
        initclock = getInitClock(data[nT / 2:nT / 2 + 100][:, :, RSstations + otherstations], freqs)  # only on a few timestamps
        logging.info('initial clocks: ' + str(initclock[1]))
        # init CS clocks to 0
        # logging.info("data before init clock" + str(data[nT/2,:,-1]))
        data[:, :, RSstations + otherstations] = data[:, :, RSstations + otherstations] - freqs[np.newaxis, :, np.newaxis, np.newaxis] * initclock[1][np.newaxis, np.newaxis, :] * 1e-9 * 2 * np.pi
        # logging.info("clock correction" + str(np.remainder(freqs*initclock[1][-1]*-1e-9*2*np.pi+np.pi,2*np.pi)-np.pi))
        # logging.info("data after init clock" + str(np.remainder(data[nT/2,:,-1]+np.pi,2*np.pi)-np.pi))
    offset = np.zeros((nSt, npol), dtype=np.float32)
    # initialize arrays
    clock = np.zeros((nT, nSt, npol), dtype=np.float32)
    tec = np.zeros((nT, nSt, npol), dtype=np.float32)
    # better not to use fitoffset
    for pol in xrange(npol):
        # get a good guesss without offset
        # logging.info("sending masked data "+str(data[:,:,:,pol].count()))
        initialchi2cut = chi2cut  # user defined
        if removePhaseWraps:
            initialchi2cut = 30000.  # this number is quite arbitrary
        (tecarray, clockarray, residualarray) = getClockTECFit(
            np.ma.copy(data[:, :, :, pol]),
            freqs,
            stations,
            initSol=initSol,
            returnResiduals=True,
            chi2cut=initialchi2cut,
            )
        if removePhaseWraps:
            # correctfrist times only,try to make init correct ?
            (offset[:, pol], wraps, steps) = correctWraps(tecarray, residualarray, freqs, station_positions)
        logging.info('residual iter 1, pol %d: ' % pol + str(residualarray[0, 0]))
        logging.info('tec iter 1, pol %d: ' % pol + str(tecarray[0]))
        logging.info('clock iter 1, pol %d: ' % pol + str(clockarray[0]))
        if removePhaseWraps:
            logging.info('wraps: ' + str(wraps))
        logging.info('offsets: ' + str(offset[:, pol]))

        data[:, :, :, pol] += offset[:, pol][np.newaxis, np.newaxis]

        # remove fitoffset
        if removePhaseWraps:
            initsol = np.zeros((nSt, 2), dtype=np.float32)
            initsol[:, 0] = tecarray[0, :] + wraps * steps[0]
            initsol[:, 1] = clockarray[0, :] + wraps * steps[1]
            logging.info('initsol tec, pol %d: ' % pol + str(initsol[:, 0]))
            logging.info('initsol clock, pol %d: ' % pol + str(initsol[:, 1]))
            tecarray = 0
            clockarray = 0
            residualarray = 0
            # is it needed to redo the fitting? this is the time bottleneck
            (tec[:, :, pol], clock[:, :, pol]) = getClockTECFit(
                np.ma.copy(data[:, :, :, pol]),
                freqs,
                stations,
                initSol=initsol,
                returnResiduals=False,
                chi2cut=chi2cut,
                )
        else:
            tec[:, :, pol] = tecarray[:, :]
            clock[:, :, pol] = clockarray[:, :]
        logging.info('tec iter 2, pol %d: ' % pol + str(tec[0, :, pol]))
        logging.info('clock iter 2, pol %d: ' % pol + str(clock[0, :, pol]))
    if not 'LBA' in stations[0] and len(initSol) < 1:
        clock[:, RSstations + otherstations] += initclock[1][np.newaxis, :, :]
    return (clock, tec, offset, stations)
