"""
This file contains tests for the losoto operation flagextend.
"""

from losoto.operations import flagextend


def test_flagextend(soltab):
    """
    Test the Losoto operation flagextend
    """
    assert (
        flagextend.run(soltab, axesToExt, size, percent=50.0, maxCycles=3, ncpu=0) == 0
    )
    pass
