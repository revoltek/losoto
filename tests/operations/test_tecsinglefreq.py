"""
This file contains tests for the losoto operation tecsinglefreq.
"""

from losoto.operations import tecsinglefreq


def test_tecsinglefreq(soltab):
    """
    Test the Losoto operation tecsinglefreq
    """
    assert tecsinglefreq.run(soltab, soltabOut="tec000", refAnt="") == 0
    pass
