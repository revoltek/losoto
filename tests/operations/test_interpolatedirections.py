"""
This file contains tests for the losoto operation interpolatedirections.
"""

from losoto.operations import interpolatedirections


def test_interpolatedirections(soltab):
    """
    Test the Losoto operation interpolatedirections
    """
    assert (
        interpolatedirections.run(
            soltab, interp_dirs, soltabOut=None, prefix="interp_", ncpu=0
        )
        == 0
    )
