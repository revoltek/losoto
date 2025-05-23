"""
This file contains tests for the losoto operation norm.
"""

from losoto.operations import norm


def test_norm(soltab):
    """
    Test the Losoto operation norm
    """
    assert norm.run(soltab, axesToNorm, normVal=1.0, log=False) == 0
