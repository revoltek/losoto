"""
This file contains tests for the losoto operation duplicate.
"""

from losoto.operations import duplicate


def test_duplicate(soltab):
    """
    Test the Losoto operation duplicate
    """
    assert duplicate.run(soltab, soltabOut="", overwrite=False) == 0
