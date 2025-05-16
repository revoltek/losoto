"""
This file contains tests for the losoto operation replicateonaxis.
"""

from losoto.operations import replicateonaxis


def test_replicateonaxis(soltab):
    """
    Test the Losoto operation replicateonaxis
    """
    assert replicateonaxis.run(soltab, axisReplicate, fromCell, updateWeights=True) == 0
