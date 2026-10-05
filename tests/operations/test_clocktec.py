"""
This file contains tests for the losoto operation clocktec.
"""

import pytest

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
            refAnt="",
            numParams=0,
            combinePol=False,
            circular=False,
            weightMode="kappa",
            initMode="linear",
            initWithPrevious=True,
            removeSteps=False,
            spatialFit=False,
            allowPiWraps=False,
            smoothSolutions=False,
            lofar2Mode=True,
            finalPh0Fit=0,
            logFile="",
        )
        == 0
    )
    solset = soltab.getSolset()
    assert solset.getSoltab("tec000").getType() == "tec"
    assert solset.getSoltab("clock000").getType() == "clock"
    assert solset.getSoltab("phase_offset000").getType() == "phase"


def test_clocktec_old_parameters(soltab):
    """
    Parameters of the old implementation must no longer be accepted
    """
    with pytest.raises(TypeError):
        clocktec.run(soltab, flagBadChannels=True)
