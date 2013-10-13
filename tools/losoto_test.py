#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This tool tests the functionalities of H5parm libraries

# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)"

import sys, os, time
import numpy as np
import logging
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm, solFetcher, solWriter

if os.path.isfile('test.h5'): os.system('rm test.h5')

# general h5parm library
logging.info("Create a new H5parm")
H5 = h5parm('test.h5', readonly=False)
logging.info("Close H5parm")
del H5
logging.info("Open in read-only mode")
H5 = h5parm('test.h5', readonly=True)
del H5
print "###########################################"

# solsets
H5 = h5parm('test.h5', readonly=False)
logging.info("Create solsets (using same names)")
H5.makeSolset('ssTest')
H5.makeSolset('ssTest')
H5.makeSolset()
logging.info('Get a solset object')
ss=H5.getSolset('ssTest')
logging.info('Get ants')
ant=H5.getAnt(ss)
logging.info('Get sources')
sou=H5.getSou(ss)
logging.info('Get all solsets:')
print H5.getSolsets()
print "###########################################"

# soltabs
logging.info("Create soltabs (using same names)")
axesVals = [['a','b','c','d'], np.arange(100), np.arange(1000)]
vals = np.arange(4*100*1000).reshape(4,100,1000)
H5.makeSoltab(ss, 'amplitude', 'stTest', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
H5.makeSoltab(ss, 'amplitude', 'stTest', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
H5.makeSoltab(ss, 'amplitude', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
logging.info('Get a soltab object')
st=H5.getSoltab(ss,'stTest')
logging.info('Get all soltabs:')
print H5.getSoltabs(ss)

print "###########################################"
logging.info('printInfo()')
print H5.printInfo()

print "###########################################"
logging.info('solFetcher/solWriter')
Hsf = solFetcher(st)
Hsw = solWriter(st)
logging.info('Get Axes Names')
print Hsf.getAxesNames()
logging.info('Get Axes1 Len (exp 4)')
print Hsf.getAxisLen('axis1')
logging.info('Get solution Type (exp: amplitude)')
print Hsf.getType()
logging.info('Get axisValues (exp: a,b,c,d)')
print Hsf.getAxisValues('axis1')
logging.info('Set new axes values')
print Hsw.setAxisValues('axis1',['e','f','g','h'])
logging.info('Get axisValues (exp: e,f,g,h)')
print Hsf.getAxisValues('axis1')
logging.info('Set a selection using single/multiple vals and append (exp: 1x1x2)')
Hsf.setSelection(axis1='e', axis2=1, axis3=[1,10])
v,a = Hsf.getValues()
print v.shape
print a
logging.info('Set a selection using min max (exp: 4x10x10)')
Hsf.setSelection(axis2={'min':10,'max':19}, axis3={'min':990, 'max':1e6})
v,a = Hsf.getValues()
print v.shape
print a
logging.info('Writing back with selction')
Hsw.setSelection(axis1='e', axis2={'min':10,'max':19}, axis3={'min':990, 'max':1e6})
Hsw.setValues(v[0])
logging.info('Set/Get history')
Hsw.addHistory('history is working.')
print Hsw.getHistory()
logging.info('Get Vaues Iter (exp: 40 and 10)')
i=0
for matrix, coord in Hsf.getValuesIter(returnAxes=['axis3']):
    i += 1
print "Iterations:", i
i=0
for matrix, coord in Hsf.getValuesIter(returnAxes=['axis2','axis3']):
    print coord
    i += 1
print "Iterations:", i


# close the H5parm file
del Hsf
del Hsw
del H5
