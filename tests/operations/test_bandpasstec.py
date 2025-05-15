"""
This file contains tests for the losoto operation bandpasstec.
"""

from losoto.operations import bandpasstec


def test_bandpasstec(soltab):
    """
    Test the Losoto operation bandpasstec
    """
    assert (
        bandpasstec.run(
            soltab, soltabOutTEC="tec000", soltabOutBP="phase000", refAnt=""
        )
        == 0
    )
