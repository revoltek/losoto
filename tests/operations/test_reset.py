"""
This file contains tests for the losoto operation reset.
"""

from losoto.operations import reset


def test_reset(soltab):
    """
    Test the Losoto operation reset
    """
    assert reset.run(soltab, dataVal=-999.0) == 0
