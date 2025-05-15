"""
This file contains tests for the losoto operation flagstation.
"""

from losoto.operations import flagstation


def test_flagstation(soltab):
    """
    Test the Losoto operation flagstation
    """
    assert (
        flagstation.run(
            soltab,
            mode,
            minFlaggedFraction=0.0,
            maxFlaggedFraction=0.5,
            nSigma=5.0,
            maxStddev=None,
            ampRange=None,
            thresholdBaddata=0.5,
            telescope="lofar",
            skipInternational=False,
            skipAnts=[],
            refAnt="",
            soltabExport="",
            ncpu=0,
        )
        == 0
    )
    pass
