#!/usr/bin/python
# -*- coding: utf-8 -*-

# Cyril Tasse
# Francesco de Gasperin
# Reinout van Weeren

import losoto.operations
import losoto._version
import lofar.parameterset
import sys, os

# Options
import optparse
opt = optparse.OptionParser(usage='%prog parset [default: losoto.parset]', version='%prog '+losoto._version.__version__)
opt.add_option('-v', '--verbose', help='GO VERBOSE! (default=False)', action='store_true', default=False)
opt.add_option('-d', '--gds', help='Gds file to construct the parmdb and the hdf5', type='string', default=None)
opt.add_option('-g', '--globaldb', help='Globaldb in hdf5 format', type='string', default=None)
(options, args) = opt.parse_args()

verb = options.verbose
gdsfile = options.gds
globaldbFile = options.globaldb

try: parsetFile = args[0]
except: parsetFile = 'losoto.parset'

if not os.path.isfile(parsetFile):
    print "Error: missing parset, I don't know what to do :'("
    sys.exit(1)

# from ~vdtol/Expion-2011-05-03/src
parset = lofar.parameterset.parameterset( parsetFile )
steps = parset.getStringVector( "LoSoTo.Steps", [] )
gainType = parset.getString( "LoSoTo.SolType", 'Gain' )
phasorsEnable = parset.getBool( "LoSoTo.Phasors.Enable", False )
antenna = parset.getStringVector( "LoSoTo.Antenna", [] )
polarizations = parset.getStringVector("LoSoTo.Polarizations", ['XX', 'YY'])

# possible operations, linked to relative function
operations = { "RESET": losoto.operations.reset,
               "CLOCKTEC": losoto.operations.clocktec,
               "FLAG": losoto.operations.flag,
               "SMOOTH": losoto.operations.smooth,
               "INTERP": losoto.operations.interp,
#               "WRITE": losoto.operations.write,
               "PLOT": losoto.operations.plot,
               "APPLY": losoto.operations.apply }

# open the HDF5
H=ClassMakeHDF5.ClassMakeHDF5()
H.load_globaldb(globaldbFile)

for step in steps:
   operation = parset.getString( '.'.join( [ "LoSoTo.Steps", step, "Operation" ] ) )
   print "Starting \'" + step + "\' step (operation: " + operation + ")."
   H = operations[ operation ] ( step, parset, H )

print "Done."
