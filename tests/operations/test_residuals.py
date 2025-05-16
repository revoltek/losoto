"""
This file contains tests for the losoto operation residuals.
"""

from losoto.operations import residuals


def test_residuals(soltab):
    """
    Test the Losoto operation residuals
    """
    assert residuals.run(soltab, soltabsToSub, ratio=False) == 0
