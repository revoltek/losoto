"""
This file contains tests for the losoto operation flag.
"""

from losoto.operations import flag


def test_flag(soltab):
    """
    Test the Losoto operation flag
    """
    assert (
        flag.run(
            soltab,
            axesToFlag,
            order,
            maxCycles=5,
            maxRms=5.0,
            maxRmsNoise=0.0,
            fixRms=0.0,
            fixRmsNoise=0.0,
            windowNoise=11,
            replace=False,
            preflagzeros=False,
            mode="smooth",
            refAnt="",
            ncpu=0,
        )
        == 0
    )
