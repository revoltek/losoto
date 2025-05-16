"""
This file contains tests for the losoto operation example.
"""

from losoto.operations import example


def test_example(soltab):
    """
    Test the Losoto operation example
    """
    assert example.run(soltab, opt1, opt2=[1.0, 2.0, 3.0], opt3=0) == 0
