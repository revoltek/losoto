"""
This file contains tests for the losoto operation tecjump.
"""

from losoto.operations import tecjump


def test_tecjump(soltab):
    """
    Test the Losoto operation tecjump
    """
    assert tecjump.run(soltab, refAnt="", soltabError="", ncpu=0) == 0
