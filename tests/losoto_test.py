#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

# This tool tests the functionalities of H5parm libraries

# Authors:
# Francesco de Gasperin
_author = "Francesco de Gasperin (fdg@hs.uni-hamurg.de)"

import sys, os, time
import numpy as np
import logging
import losoto._version
import losoto._logging
from losoto.h5parm import h5parm

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
print("###########################################")
logging.info('### Solset')
H5 = h5parm('test.h5', readonly=False)
logging.info("Create solset")
H5.makeSolset('ssTest')
logging.info("Create solsets (using same name)")
H5.makeSolset('ssTest')
logging.info("Create solsets (using default name)")
ssdel = H5.makeSolset()
logging.info("Delete solset")
ssdel.delete()
logging.info('Get all solsets:')
print(H5.getSolsetNames())
logging.info('Get a solset object')
ss=H5.getSolset('ssTest')
logging.info('Get ants')
ant=ss.getAnt()
logging.info('Get sources')
sou=ss.getSou()

# soltabs
print("###########################################")
logging.info('### Soltab')
axesVals = [['a','b','c','d'], np.arange(10, dtype=float), np.arange(100, dtype=float)]
vals = np.arange(4*10*100).reshape(4,10,100)
logging.info("Create soltab")
ss.makeSoltab('amplitude', 'stTest', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
logging.info("Create soltab (using same name)")
ss.makeSoltab('amplitude', 'stTest', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
logging.info("Create soltab (using default name)")
stdel = ss.makeSoltab('amplitude', axesNames=['axis1','axis2','axis3'], axesVals=axesVals, vals=vals, weights=vals)
logging.info('Del Soltab')
stdel.delete()
logging.info('Get a soltab object')
st=ss.getSoltab('stTest')
logging.info('Get all soltabs:')
print(ss.getSoltabNames())

print("###########################################")
logging.info('### Soltab - R/W')
logging.info('Get solution Type (exp: amplitude)')
print(st.getType())
logging.info('Get Axes Names')
print(st.getAxesNames())
logging.info('Get Axis1 Len (exp: 4)')
print(st.getAxisLen('axis1'))
logging.info('Get Axis1 Type (exp: str)')
print(st.getAxisType('axis1'))
logging.info('Get Axis2 Type (exp: float)')
print(st.getAxisType('axis2'))
logging.info('Get Axis1 Values (exp: a,b,c,d)')
print(st.getAxisValues('axis1'))
logging.info('Set new axes values')
st.setAxisValues('axis1',['e','f','g','h'])
logging.info('Get new Axis1 Values (exp: e,f,g,h)')
print(st.getAxisValues('axis1'))

print("###########################################")
logging.info('### Soltab - selection')
logging.info('Set a selection using single/multiple vals and append (exp: 3x1x2)')
st.setSelection(axis1=['f','e','h'], axis2=1., axis3=[1.,10.])
v,a = st.getValues()
print(v.shape)
print(v)
print(a)
logging.info('Writing back with selction')
st.setSelection(axis1=['f','e','h'], axis2=1., axis3=[1.,10.])
st.setValues(v)

logging.info('Set a selection using min max (exp: 2x4x10)')
st.setSelection(axis1=['e','h'], axis2={'min':2,'max':5}, axis3={'min':90, 'max':1e6})
v,a = st.getValues()
print(v.shape)
print(a)
logging.info('Writing back with selction')
st.setSelection(axis1=['e','h'], axis2={'min':2,'max':5}, axis3={'min':90, 'max':1e6})
st.setValues(v)

logging.info('Get Vaues Iter (exp: 10)')
i=0
for matrix, coord, sel in st.getValuesIter(returnAxes=['axis3']):
    i += 1
print(matrix.shape)
print("Iterations:", i, "(expected: 2x4=8)")
logging.info('Get Vaues Iter (exp: 4x10)')
i=0
for matrix, coord, sel in st.getValuesIter(returnAxes=['axis2','axis3']):
    i += 1
print(matrix.shape)
print("Iterations:", i, "(expected: 2)")


print("###########################################")
logging.info('### Soltab - History and info')
logging.info('Set a selection using single/multiple vals and append (exp: 3x1x2)')
logging.info('Set/Get history')
st.addHistory('History is working.')
print(st.getHistory())

logging.info('printInfo()')
print(H5.printInfo())

# close the H5parm file
del st
del ss
del H5
os.system('rm test.h5')
logging.info('Done.')
