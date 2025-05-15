"""
This file contains tests for the losoto operation deleteaxis.
"""

from losoto.operations import deleteaxis


def test_deleteaxis(soltab):
    """
    Test the Losoto operation deleteaxis
    """
    assert deleteaxis.run(soltab, axisDelete, fromCell) == 0
    pass
