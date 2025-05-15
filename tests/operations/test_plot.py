"""
This file contains tests for the losoto operation plot.
"""

from losoto.operations import plot


def test_plot(soltab):
    """
    Test the Losoto operation plot
    """
    assert (
        plot.run(
            soltab,
            axesInPlot,
            axisInTable="",
            axisInCol="",
            axisDiff="",
            NColFig=0,
            figSize=[0, 0],
            markerSize=2,
            minmax=[0, 0],
            log="",
            plotFlag=False,
            doUnwrap=False,
            refAnt="",
            refDir="",
            soltabsToAdd="",
            makeAntPlot=False,
            makeMovie=False,
            prefix="",
            ncpu=0,
        )
        == 0
    )
    pass
