"""
This file contains tests for the losoto operation prefactor_bandpass.
"""

from losoto.operations import prefactor_bandpass


def test_prefactor_bandpass(soltab):
    """
    Test the Losoto operation prefactor_bandpass
    """
    assert (
        prefactor_bandpass.run(
            soltab,
            chanWidth="",
            outSoltabName="bandpass",
            BadSBList="",
            interpolate=True,
            removeTimeAxis=True,
            autoFlag=False,
            nSigma=5.0,
            maxFlaggedFraction=0.5,
            maxStddev=0.01,
            ncpu=0,
        )
        == 0
    )
    pass
