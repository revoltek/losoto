"""
This file contains tests for the losoto operation screenvalues.
"""

from losoto.operations import screenvalues


def test_screenvalues(soltab):
    """
    Test the Losoto operation screenvalues
    """
    assert screenvalues.run(soltab1, source_dict, outsoltab, soltab2=None, ncpu=0) == 0
    pass
