"""
This file contains tests for the losoto operation plotscreen.
"""

from losoto.operations import plotscreen


def test_plotscreen(soltab):
    """
    Test the Losoto operation plotscreen
    """
    assert (
        plotscreen.run(
            soltab,
            resSoltab="",
            minZ=-3.2,
            maxZ=3.2,
            prefix="",
            remove_gradient=False,
            show_source_names=False,
            ncpu=0,
        )
        == 0
    )
