"""
This file contains tests for the losoto operation interpolate.
"""

from losoto.operations import interpolate


def test_interpolate(soltab):
    """
    Test the Losoto operation interpolate
    """
    assert (
        interpolate.run(
            soltab,
            outsoltab,
            axisToRegrid,
            newdelta,
            delta="",
            maxFlaggedWidth=0,
            log=False,
        )
        == 0
    )
