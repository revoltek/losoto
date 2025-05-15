"""
This file contains tests for the losoto operation frjump.
"""

from losoto.operations import frjump


def test_frjump(soltab):
    """
    Test the Losoto operation frjump
    """
    assert frjump.run(soltab, soltabOut, clipping, soltabPhase, frequencies) == 0
    pass
