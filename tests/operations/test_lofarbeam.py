"""
This file contains tests for the losoto operation lofarbeam.
"""

from losoto.operations import lofarbeam


def test_lofarbeam(soltab):
    """
    Test the Losoto operation lofarbeam
    """
    assert (
        lofarbeam.run(
            soltab, ms, useElementResponse=True, useArrayFactor=True, allChan=False
        )
        == 0
    )
    pass
