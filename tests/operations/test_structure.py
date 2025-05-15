"""
This file contains tests for the losoto operation structure.
"""

from losoto.operations import structure


def test_structure(soltab):
    """
    Test the Losoto operation structure
    """
    assert structure.run(soltab, doUnwrap=False, refAnt="", plotName="", ndiv=1) == 0
    pass
