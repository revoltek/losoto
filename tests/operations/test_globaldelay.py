"""
This file contains tests for the losoto operation globaldelay.
"""

from losoto.operations import globaldelay


def test_globaldelay(soltab):
    """
    Test the Losoto operation globaldelay
    """
    assert globaldelay.run(soltab, soltabOut="tec000", refAnt="") == 0
