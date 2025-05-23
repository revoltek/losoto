"""
This module contains unit tests for the h5parm module.

TODO: Here's the complete list of classes/functions in h5parm.py.
      Some of these, we may want to test too.

    def openSoltab(h5parmFile, solsetName=None, soltabName=None, address=None, readonly=True):
    class h5parm( object ):
        def __init__(self, h5parmFile, readonly=True, complevel=0, complib='zlib'):
        def open(self):
        def close(self):
        def __enter__(self):
        def __exit__(self, *exc):
        def __str__(self):
        def makeSolset(self, solsetName=None, addTables=True):
        def getSolsets(self):
        def getSolsetNames(self):
        def getSolset(self, solset):
        def _firstAvailSolsetName(self):
        def printInfo(self, filter=None, verbose=False):
            def grouper(n, iterable, fillvalue=' '):
            def wrap(text, width=80):
    class Solset( object ):
        def __init__(self, solset):
        def close(self):
        def delete(self):
        def rename(self, newname, overwrite=False):
        def makeSoltab(self, soltype=None, soltabName=None,
        def _fisrtAvailSoltabName(self, soltype):
        def getSoltabs(self, useCache=False, sel={}):
        def getSoltabNames(self):
        def getSoltab(self, soltab, useCache=False, sel={}):
        def getAnt(self):
        def getSou(self):
        def getAntDist(self, ant=None, ant_subset=[]):
    class Soltab( object ):
        def __init__(self, soltab, useCache = False, args = {}):
        def delete(self):
        def rename(self, newname, overwrite=False):
        def setCache(self, val, weight):
        def getSolset(self):
        def getAddress(self):
        def clearSelection(self):
        def setSelection(self, update=False, **args):
        def getType(self):
        def getAxesNames(self):
        def getAxisLen(self, axis, ignoreSelection=False):
        def getAxisType(self, axis):
        def getAxisValues(self, axis, ignoreSelection=False):
        def setAxisValues(self, axis, vals):
        def setValues(self, vals, selection = None, weight = False):
        def flush(self):
        def __getattr__(self, axis):
        def _applyAdvSelection(self, data, selection, ignoreSelection=False):
        def _getFullyFlaggedAnts(self):
        def getValues(self, retAxesVals=True, weight=False, refAnt=None, refDir=None):
        def getValuesIter(self, returnAxes=[], weight=False, refAnt=None, refDir=None):
            def g():
        def autoRefAnt(self):
        def addHistory(self, entry, date = True):
        def getHistory(self):
"""

import numpy
import os
import pytest
import tempfile
from losoto.h5parm import h5parm


@pytest.fixture
def h5parm_file(tmp_path):
    """Fixture to create a temporary h5parm file."""
    fname = tempfile.mktemp(suffix=".h5")
    yield fname
    os.remove(fname)


def test_h5parm(h5parm_file):
    """
    Test the h5parm module.
    TODO: This is a very basic test. We should add more tests to cover
    different scenarios and functionalities of the h5parm class.
    """

    with h5parm(h5parm_file, readonly=False) as h5:
        solset = h5.makeSolset("sol000")
        dirvals = ["CasA", "VirA"]
        timevals = numpy.arange(0, 5)
        freqvals = numpy.arange(0, 3)
        vals = numpy.ones((2, 3, 5))
        soltab = solset.makeSoltab(
            soltype="phase",
            soltabName="phase000",
            axesNames=["dir", "freq", "time"],
            axesVals=[dirvals, freqvals, timevals],
            vals=vals,
            weights=vals,
        )

    with h5parm(h5parm_file, readonly=True) as h5:
        solset = h5.getSolset("sol000")
        soltab = solset.getSoltab("phase000")
        vals = soltab.getValues(retAxesVals=False)
        axesnames = soltab.getAxesNames()

        assert vals.shape == (2, 3, 5)
        assert axesnames == ["dir", "freq", "time"]


def test_solset():
    """
    Test the Solset class.
    """
    pass


def test_soltab():
    """
    Test the Soltab class.
    """
    pass
