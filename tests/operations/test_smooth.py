"""
This file contains tests for the losoto operation smooth.
"""

from losoto.operations import smooth


def test_smooth(soltab):
    """
    Test the Losoto operation smooth
    """
    assert (
        smooth.run(
            soltab,
            axesToSmooth,
            size=[],
            mode="runningmedian",
            degree=1,
            replace=False,
            log=False,
            refAnt="",
        )
        == 0
    )
