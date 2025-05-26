"""
This file contains tests for the losoto operation tec.
"""

from losoto.operations import tec


def test_tec(soltab):
    """
    Test the Losoto operation tec
    """
    assert (
        tec.run(soltab, soltabOut, refAnt, maxResidualFlag, maxResidualProp, ncpu) == 0
    )
