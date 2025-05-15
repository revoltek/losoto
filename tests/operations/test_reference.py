"""
This file contains tests for the losoto operation reference.
"""

from losoto.operations import reference


def test_reference(soltab):
    """
    Test the Losoto operation reference
    """
    assert reference.run(soltab, refAnt="", refDir="") == 0
    pass
