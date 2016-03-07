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
H5.close()
logging.info("Open in read-only mode")
H5 = h5parm('test.h5', readonly=True)
H5.close()

# solsets
print "###########################################"
logging.info('### solSets')
H5 = h5parm('test.h5', readonly=False)
logging.info("Create solset")
H5.makeSolset('ssTest')
logging.info("Create solsets (using same name)")
H5.makeSolset('ssTest')
logging.info("Create solsets (using default name)")
H5.makeSolset()
logging.info('Get a solset object')
ss=H5.getSolset('ssTest')
logging.info('Get ants')
ant=H5.getAnt(ss)
logging.info('Get sources')
sou=H5.getSou(ss)
logging.info('Get all solsets:')
print H5.getSolsets()

# soltabs
print "###########################################"
logging.info('### solTabs')
axesVals = [['a','b','c','d'], np.arange(10), np.arange(100)]
vals = np.arange(4*10*100).reshape(4,10,100)
logging.info("Create soltab")
H5.makeSoltab(ss, 'amplitude', 'stTest', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
logging.info("Create soltab (using same name)")
H5.makeSoltab(ss, 'amplitude', 'stTest', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
logging.info("Create soltab (using default name)")
H5.makeSoltab(ss, 'amplitude', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
logging.info('Get a soltab object')
st=H5.getSoltab(ss,'stTest')
logging.info('Get all soltabs:')
print H5.getSoltabs(ss)

print "###########################################"
logging.info('### solFetcher/solWriter - General')
Hsf = solFetcher(st)
Hsw = solWriter(st)
logging.info('Get solution Type (exp: amplitude)')
print Hsf.getType()
logging.info('Get Axes Names')
print Hsf.getAxesNames()
logging.info('Get Axis1 Len (exp: 4)')
print Hsf.getAxisLen('axis1')
logging.info('Get Axis1 Type (exp: str)')
print Hsf.getAxisType('axis1')
logging.info('Get Axis2 Type (exp: float)')
print Hsf.getAxisType('axis2')
logging.info('Get Axis1 Values (exp: a,b,c,d)')
print Hsf.getAxisValues('axis1')
logging.info('Set new axes values')
Hsw.setAxisValues('axis1',['e','f','g','h'])
logging.info('Get new Axis1 Values (exp: e,f,g,h)')
print Hsf.getAxisValues('axis1')

print "###########################################"
logging.info('### solFetcher/solWriter - selection')
logging.info('Set a selection using single/multiple vals and append (exp: 3x1x2)')
Hsf.setSelection(axis1=['f','e','h'], axis2=1, axis3=[1,10])
v,a = Hsf.getValues()
print v.shape
print v
print a
logging.info('Writing back with selction')
Hsw.setSelection(axis1=['f','e','h'], axis2=1, axis3=[1,10])
Hsw.setValues(v)

logging.info('Set a selection using min max (exp: 4x4x10)')
Hsf.setSelection(axis1='e', axis2={'min':2,'max':5}, axis3={'min':90, 'max':1e6})
v,a = Hsf.getValues()
print v.shape
print a
logging.info('Writing back with selction')
Hsw.setSelection(axis1='e', axis2={'min':2,'max':5}, axis3={'min':90, 'max':1e6})
Hsw.setValues(v)

logging.info('Get Vaues Iter (exp: 40 and 10)')
i=0
for matrix, coord, sel in Hsf.getValuesIter(returnAxes=['axis3']):
    print matrix.shape
    i += 1
print "Iterations:", i
i=0
for matrix, coord, sel in Hsf.getValuesIter(returnAxes=['axis2','axis3']):
    print matrix.shape
    print coord
    i += 1
print "Iterations:", i


print "###########################################"
logging.info('### solHandler - History and info')
logging.info('Set a selection using single/multiple vals and append (exp: 3x1x2)')
logging.info('Set/Get history')
Hsw.addHistory('History is working.')
print Hsw.getHistory()

logging.info('printInfo()')
print H5.printInfo()

# close the H5parm file
del Hsf
del Hsw
del H5
os.system('rm test.h5')
logging.info('Done.')
