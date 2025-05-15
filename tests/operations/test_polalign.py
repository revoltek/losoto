"""
This file contains tests for the losoto operation polalign.
"""

from losoto.operations import polalign


def test_polalign(soltab):
    """
    Test the Losoto operation polalign
    """
    assert polalign.run(soltab, soltabOut="phasediff", minFreq=0, refAnt="") == 0
