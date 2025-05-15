"""
This file contains tests for the Losoto operations.
Each function tests a specific operation in the Losoto package.
"""

import numpy
import os
import pytest
import tempfile

from losoto import operations as op
from losoto.h5parm import h5parm
from losoto._logging import logger as logging


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
    solset = h5.makeSolset("sol000")
    phases = numpy.array([2, 3, 5])
    weights = numpy.ones(len(phases))
    freqs = numpy.arange(len(phases))
    return solset.makeSoltab(
        soltype="amplitude",
        soltabName="phase000",
        axesNames=["phase"],
        axesVals=[freqs],
        vals=phases,
        weights=weights,
    )  # , axes=['time', 'frequency', 'antenna'])
