"""
This file contains tests for the losoto operation stationscreen.
"""

from losoto.operations import stationscreen


def test_stationscreen(soltab):
    """
    Test the Losoto operation stationscreen
    """
    assert (
        stationscreen.run(
            soltab,
            outsoltab,
            order=12,
            beta=5.0 / 3.0,
            niter=2,
            nsigma=5.0,
            refAnt=-1,
            scale_order=True,
            scale_dist=None,
            min_order=5,
            adjust_order=True,
            ncpu=0,
        )
        == 0
    )
