"""
This file contains tests for the losoto operation clip.
"""

from losoto.operations import clip


def test_clip(soltab):
    """
    Test the Losoto operation clip
    """
    assert (
        clip.run(soltab, axesToClip=None, clipLevel=5.0, log=False, mode="median") == 0
    )
    pass
