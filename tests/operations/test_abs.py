"""
This file contains tests for the losoto operation abs.
"""

from losoto.operations import abs


def test_abs(soltab):
    """
    Test the Losoto operation abs
    """
    assert abs.run(soltab) == 0
