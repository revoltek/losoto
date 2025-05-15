"""
This file contains tests for the losoto operation interpolate.
"""

from losoto.operations import interpolate


def test_interpolate(soltab):
    """
    Test the Losoto operation interpolate
    """
    assert (
        interpolate.run(ssoltab, interp_dirs, soltabOut=None, prefix="interp_", ncpu=0)
        == 0
    )
    pass
