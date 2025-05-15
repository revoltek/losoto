"""
This file contains tests for the losoto operation directionscreen.
"""

from losoto.operations import directionscreen


def test_directionscreen(soltab):
    """
    Test the Losoto operation directionscreen
    """
    assert (
        directionscreen.run(
            soltab,
            outSoltab="tecscreen",
            height=200.0e3,
            order=12,
            beta=5.0 / 3.0,
            ncpu=0,
        )
        == 0
    )
    pass
