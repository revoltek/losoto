"""
This file contains tests for the Losoto operations.
Each function tests a specific operation in the Losoto package.
"""
import os
import pytest
import tempfile

from losoto import operations as op
from losoto.h5parm import h5parm

@pytest.fixture
def h5():
    """Fixture to create a temporary h5parm file."""
    fname = tempfile.mktemp(suffix=".h5")
    h5p = h5parm(fname, readonly=False)
    yield h5p
    h5p.close()
    os.remove(fname)

@pytest.fixture
def soltab(h5):
    """
    Fixture to create a sample solution table for testing.
    This should be replaced with actual code to create a solution table.
    """
    solset = h5.makeSolset('sol000')
    return solset.makeSoltab()

def test_abs(soltab):
    """
    Test the Losoto operation abs
    """
    op.abs.run(soltab)
    pass

def test_bandpasstec(soltab):
    """
    Test the Losoto operation bandpasstec
    """
    op.bandpasstec.run(soltab, soltabOutTEC='tec000', soltabOutBP='phase000', refAnt='')
    pass

def test_clip(soltab):
    """
    Test the Losoto operation clip
    """
    op.clip.run(soltab, axesToClip=None, clipLevel=5., log=False, mode='median')
    pass

def test_clocktec(soltab):
    """
    Test the Losoto operation clocktec
    """
    op.clocktec.run(soltab, tecsoltabOut='tec000', clocksoltabOut='clock000', offsetsoltabOut='phase_offset000', tec3rdsoltabOut='tec3rd000', flagBadChannels=True, flagCut=5., chi2cut=3000., combinePol=False, removePhaseWraps=True, fit3rdorder=False, circular=False, reverse=False, invertOffset=False, nproc=10)
    pass

def test_deleteaxis(soltab):
    """
    Test the Losoto operation deleteaxis
    """
    op.deleteaxis.run(soltab, axisDelete, fromCell)
    pass

def test_directionscreen(soltab):
    """
    Test the Losoto operation directionscreen
    """
    op.directionscreen.run(soltab, outSoltab='tecscreen', height=200.0e3, order=12,
    beta=5.0/3.0, ncpu=0)
    pass

def test_duplicate(soltab):
    """
    Test the Losoto operation duplicate
    """
    op.duplicate.run(soltab, soltabOut='', overwrite=False)
    pass

def test_example(soltab):
    """
    Test the Losoto operation example
    """
    op.example.run(soltab, opt1, opt2 = [1., 2., 3.], opt3 = 0)
    pass

def test_faraday(soltab):
    """
    Test the Losoto operation faraday
    """
    op.faraday.run(soltab, soltabOut='rotationmeasure000', refAnt='', maxResidual=1.,ncpu=0)
    pass

def test_flagextend(soltab):
    """
    Test the Losoto operation flagextend
    """
    op.flagextend.run(soltab, axesToExt, size, percent=50., maxCycles=3, ncpu=0)
    pass

def test_flag(soltab):
    """
    Test the Losoto operation flag
    """
    op.flag.run(soltab, axesToFlag, order, maxCycles=5, maxRms=5., maxRmsNoise=0., fixRms=0., fixRmsNoise=0., windowNoise=11, replace=False, preflagzeros=False, mode='smooth', refAnt='', ncpu=0)
    pass

def test_flagstation(soltab):
    """
    Test the Losoto operation flagstation
    """
    op.flagstation.run(soltab, mode, minFlaggedFraction=0.0, maxFlaggedFraction=0.5, nSigma=5.0,
         maxStddev=None, ampRange=None, thresholdBaddata = 0.5, telescope='lofar',
         skipInternational=False, skipAnts=[], refAnt='', soltabExport='', ncpu=0)
    pass

def test_frjump(soltab):
    """
    Test the Losoto operation frjump
    """
    op.frjump.run(soltab, soltabOut,clipping,soltabPhase,frequencies)
    pass

def test_globaldelay(soltab):
    """
    Test the Losoto operation globaldelay
    """
    op.globaldelay.run(soltab, soltabOut='tec000', refAnt='')
    pass

def test_interpolatedirections(soltab):
    """
    Test the Losoto operation interpolatedirections
    """
    op.interpolatedirections.run(soltab=None)
    pass

def test_interpolate(soltab):
    """
    Test the Losoto operation interpolate
    """
    op.interpolate.run(ssoltab, interp_dirs, soltabOut=None, prefix='interp_', ncpu=0)
    pass

def test_lofarbeam(soltab):
    """
    Test the Losoto operation lofarbeam
    """
    op.lofarbeam.run(soltab, ms, useElementResponse=True, useArrayFactor=True, allChan=False)
    pass

def test_norm(soltab):
    """
    Test the Losoto operation norm
    """
    op.norm.run(soltab, axesToNorm, normVal = 1., log = False)
    pass

def test_plot(soltab):
    """
    Test the Losoto operation plot
    """
    op.plot.run(soltab, axesInPlot, axisInTable='', axisInCol='', axisDiff='', NColFig=0, figSize=[0,0], markerSize=2, minmax=[0,0], log='', \
               plotFlag=False, doUnwrap=False, refAnt='', refDir='', soltabsToAdd='', makeAntPlot=False, makeMovie=False, prefix='', ncpu=0)
    pass

def test_plotscreen(soltab):
    """
    Test the Losoto operation plotscreen
    """
    op.plotscreen.run(soltab, resSoltab='', minZ=-3.2, maxZ=3.2, prefix='', remove_gradient=False,
    show_source_names=False, ncpu=0)
    pass

def test_polalign(soltab):
    """
    Test the Losoto operation polalign
    """
    op.polalign.run(soltab, soltabOut='phasediff', minFreq=0, refAnt='')
    pass

def test_prefactor_bandpass(soltab):
    """
    Test the Losoto operation prefactor_bandpass
    """
    op.prefactor_bandpass.run(soltab, chanWidth='', outSoltabName='bandpass', BadSBList = '', interpolate=True,
        removeTimeAxis=True, autoFlag=False, nSigma=5.0, maxFlaggedFraction=0.5,
        maxStddev=0.01, ncpu=0)
    pass

def test_prefactor_XYoffset(soltab):
    """
    Test the Losoto operation prefactor_XYoffset
    """
    op.prefactor_XYoffset.run(soltab, chanWidth)
    pass

def test_reference(soltab):
    """
    Test the Losoto operation reference
    """
    op.reference.run(soltab, refAnt='', refDir='')
    pass

def test_replicateonaxis(soltab):
    """
    Test the Losoto operation replicateonaxis
    """
    op.replicateonaxis.run(soltab, axisReplicate, fromCell, updateWeights=True)
    pass

def test_reset(soltab):
    """
    Test the Losoto operation reset
    """
    op.reset.run(soltab, dataVal=-999.)
    pass

def test_residuals(soltab):
    """
    Test the Losoto operation residuals
    """
    op.residuals.run(soltab, soltabsToSub, ratio=False)
    pass

def test_reweight(soltab):
    """
    Test the Losoto operation reweight
    """
    op.reweight.run(soltab, mode='uniform', weightVal=1., nmedian=3, nstddev=251,
    soltabImport='', flagBad=False, ncpu=0)
    pass

def test_screenvalues(soltab):
    """
    Test the Losoto operation screenvalues
    """
    op.screenvalues.run(soltab1, source_dict, outsoltab, soltab2=None, ncpu=0)
    pass

def test_smooth(soltab):
    """
    Test the Losoto operation smooth
    """
    op.smooth.run(soltab, axesToSmooth, size=[], mode='runningmedian', degree=1, replace=False, log=False, refAnt='')
    pass

def test_splitleak(soltab):
    """
    Test the Losoto operation splitleak
    """
    op.splitleak.run(soltab, soltabOutG=None, soltabOutD=None)
    pass

def test_stationscreen(soltab):
    """
    Test the Losoto operation stationscreen
    """
    op.stationscreen.run(soltab, outsoltab, order=12, beta=5.0/3.0, niter=2, nsigma=5.0,
    refAnt=-1, scale_order=True, scale_dist=None, min_order=5, adjust_order=True, ncpu=0)
    pass

def test_structure(soltab):
    """
    Test the Losoto operation structure
    """
    op.structure.run(soltab, doUnwrap=False, refAnt='', plotName='', ndiv=1)
    pass

def test_tecjump(soltab):
    """
    Test the Losoto operation tecjump
    """
    op.tecjump.run(soltab, refAnt='', soltabError='', ncpu=0)
    pass

def test_tec(soltab):
    """
    Test the Losoto operation tec
    """
    op.tec.run(soltab, soltabOut, refAnt, maxResidualFlag, maxResidualProp, ncpu)
    pass

def test_tecsinglefreq(soltab):
    """
    Test the Losoto operation tecsinglefreq
    """
    op.tecsinglefreq.run(soltab, soltabOut='tec000', refAnt='')
    pass

