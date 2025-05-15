"""
This file contains tests for the losoto operation prefactor_XYoffset.
"""

from losoto.operations import prefactor_XYoffset


def test_prefactor_XYoffset(soltab):
    """
    Test the Losoto operation prefactor_XYoffset
    """
    assert prefactor_XYoffset.run(soltab, chanWidth) == 0
    pass
