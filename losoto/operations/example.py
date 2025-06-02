#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is an example operation for LoSoTo

from losoto.lib_operations import *
from losoto._logging import logger as logging

logging.debug('Loading EXAMPLE module.')

# this funct is called by losoto to set parameters and call the real run()
def _run_parser(soltab, parser, step):
    opt1 = parser.getfloat( step, 'opt1') # no default
    opt2 = parser.getarrayfloat( step, 'opt3', [1., 2., 3.])
    opt3 = parser.getint( step, 'opt2', 0 )
    
    parser.checkSpelling( step, soltab, ['op1', 'opt2', 'opt3'])
    return run(soltab, opt1, opt2, opt3)

# this funct can be called by python directly
# parameters that are non optional require the default value equal to the one defined for the parset above
def run( soltab, opt1, opt2 = [1., 2., 3.], opt3 = 0 ):
    """
    Generic unspecified step for easy expansion.

    Parameters
    ----------
    opt1 : float
        Is a mandatory parameter.

    opt2 : list of float, optional
        Is optional, by default [1.,2.,3.]

    opt2 : int, optional
        Is optional, by default 0.
    """

    # load specific libs
    import numpy as np
 
    # initial logging
    logging.info("Working on soltab: "+soltab.name)
    
    # check input
    # ...

    axisNames = soltab.getAxesNames()
    logging.info("Axis names are: "+str(axisNames))
 
    solType = soltab.getType()
    logging.info("Soltab type is: "+solType)
 
    soltab.setSelection(ant=soltab.getAxisValues('ant')[0])
    logging.info("Selection is: "+str(soltab.selection))
 
    # find axis values
    logging.info("Antennas (no selection) are: "+str(soltab.getAxisValues('ant', ignoreSelection=True)))
    logging.info("Antennas (with selection) are: "+str(soltab.getAxisValues('ant')))
    # but one can also use (selection is active here!)
    logging.info("Antennas (other method) are: "+str(soltab.ant))
    logging.info("Frequencies are: "+str(soltab.freq))
    logging.info("Directions are: "+str(soltab.dir))
    logging.info("Polarizations are: "+str(soltab.pol))
    # try to access a non-existent axis
    soltab.getAxisValues('nonexistantaxis')
 
    # now get all values given this selection
    logging.info("Get data using soltab.val")
    val = soltab.val
    logging.debug('shape of val: '+str(soltab.val.shape))
    logging.info("$ val is "+str(val[0,0,0,0,100]))
    weight = soltab.weight
    time = soltab.time
    thisTime = soltab.time[100]
 
    # another way to get the data is using the getValues()
    logging.info("Get data using getValues()")
    grid, axes = soltab.getValues()
    # axis names
    logging.info("Axes: "+str(soltab.getAxesNames()))
    # axis shape
    print(axes)
    print([soltab.getAxisLen(axis) for axis in axes]) # not ordered, is a dict!
    # data array shape (same of axis shape)
    logging.info("Shape of values: "+str(grid.shape))
    #logging.info("$ val is "+str(grid[0,0,0,0,100]))
 
    # reset selection
    soltab.setSelection()
    logging.info('Reset selection to \'\'')
    logging.info("Antennas are: "+str(soltab.ant))
    logging.info("Frequencies are: "+str(soltab.freq))
    logging.info("Directions are: "+str(soltab.dir))
    logging.info("Polarizations are: "+str(soltab.pol))
 
    # finally the getValuesIter allaws to iterate across all possible combinations of a set of axes
    logging.info('Iteration on time/freq')
    for vals, coord, selection in soltab.getValuesIter(returnAxes=['time','freq']):
        # writing back the solutions
        soltab.setValues(vals, selection)
    logging.info('Iteration on time')
    for vals, coord, selection in soltab.getValuesIter(returnAxes=['time']):
        # writing back the solutions
        soltab.setValues(vals, selection)   
    logging.info('Iteration on dir after selection to 1 dir')
    soltab.setSelection(dir='pointing') 
    for vals, coord, selection in t.getValuesIter(returnAxes=['dir']):
        # writing back the solutions
        soltab.setValues(vals, selection)
  
    return 0 # if everything went fine, otherwise 1
 
 
