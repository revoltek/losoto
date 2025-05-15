"""
This file contains tests for the losoto operation reweight.
"""

from losoto.operations import reweight


def test_reweight(soltab):
    """
    Test the Losoto operation reweight
    """
    assert (
        reweight.run(
            soltab,
            mode="uniform",
            weightVal=1.0,
            nmedian=3,
            nstddev=251,
            soltabImport="",
            flagBad=False,
            ncpu=0,
        )
        == 0
    )
    pass
