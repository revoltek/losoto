"""
This file contains tests for the losoto operation splitleak.
"""

from losoto.operations import splitleak


def test_splitleak(soltab):
    """
    Test the Losoto operation splitleak
    """
    assert splitleak.run(soltab, soltabOutG=None, soltabOutD=None) == 0
    pass
