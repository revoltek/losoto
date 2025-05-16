"""
This file contains tests for the losoto operation faraday.
"""

from losoto.operations import faraday


def test_faraday(soltab):
    """
    Test the Losoto operation faraday
    """
    assert (
        faraday.run(
            soltab, soltabOut="rotationmeasure000", refAnt="", maxResidual=1.0, ncpu=0
        )
        == 0
    )
