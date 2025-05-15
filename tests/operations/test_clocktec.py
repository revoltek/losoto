"""
This file contains tests for the losoto operation clocktec.
"""

from losoto.operations import clocktec


def test_clocktec(soltab):
    """
    Test the Losoto operation clocktec
    """
    assert (
        clocktec.run(
            soltab,
            tecsoltabOut="tec000",
            clocksoltabOut="clock000",
            offsetsoltabOut="phase_offset000",
            tec3rdsoltabOut="tec3rd000",
            flagBadChannels=True,
            flagCut=5.0,
            chi2cut=3000.0,
            combinePol=False,
            removePhaseWraps=True,
            fit3rdorder=False,
            circular=False,
            reverse=False,
            invertOffset=False,
            nproc=10,
        )
        == 0
    )
    pass
