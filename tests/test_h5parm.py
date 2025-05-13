"""
This module contains unit tests for the h5parm module.
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
