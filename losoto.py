#!/usr/bin/python
# -*- coding: utf-8 -*-

# Cyril Tasse
# Francesco de Gasperin
# Reinout van Weeren

from losoto.operations import *
import losoto._version
import lofar.parameterset
import sys, os

# Options
import optparse
opt = optparse.OptionParser(usage='%prog parset [default: losoto.parset]', version='%prog '+_version.__verions__)
#opt.add_option('-p', '--plots', help='Produces plots (default=False)', type='store_true', default=False)
opt.add_option('-v', '--verbose', help='GO VERBOSE! (default=False)', action='store_true', default=False)
opt.add_option('-d', '--gds', help='Gds file to construct the parmdb and the hdf5', type='string', default=None)
opt.add_option('-g', '--globaldb', help='Globaldb in hdf5 format', type='string', default=None)
(options, args) = opt.parse_args()

verb = options.verbose
gdsfile = options.gds
globaldbfile = options.globaldb
parsetfile = args[0]
if os.path.isfile(parsetfile):
    print "Error: missing parset, I don't know what to do :'("
    sys.exit(1)

# from ~vdtol/Expion-2011-05-03/src
parset = lofar.parameterset.parameterset( parsetfile )
steps = parset.getStringVector( "LoSoTo.Steps", [] )
DirectionalGainEnable = parset.getBool( "LoSoTo.DirectionalGain.Enable", False )
GainEnable = parset.getBool( "LoSoTo.Gain.Enable", False )
PhasorsEnable = parset.getBool( "LoSoTo.Phasors.Enable", False )
RotationEnable = parset.getBool( "LoSoTo.Rotation.Enable", False )
CommonRotationEnable = parset.getBool( "LoSoTo.CommonRotation.Enable", False )
Antenna = parset.getStringVector( "LoSoTo.Antenna", [] )
Polarizations = parset.getStringVector("LoSoTo.Polarizations", ['XX', 'YY'])

# possible operations, linked to relative function
Operations = { "SETTOZERO": operation_settozero,
               "CLOCKTEC":  operation_clocktec,
               "FLAG": operation_flag,
               "SMOOTH":  operation_smooth,
               "INTERP":  operation_interpolate,
               "WRITE":  operation_write,
               "PLOT":  operation_plot,
               "APPLY": operation_appyl }

for step in steps:
   operation = parset.getString( '.'.join( [ "LoSoTo.Steps", step, "Operation" ] ) )
   Operations[ operation ] ( step, parset )

print "Done."
