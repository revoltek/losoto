LoSoTo: LOFAR solutions tool

This library have tools to deal with LOFAR solution tables.

- fetches parmdb tables from the cluster and creates an hdf5 file
- easy access to hdf5 file
- solutions smoothing
- solutions interpolation
- set to zero solutions
- clock/tec separation


.
├── ClassMakeHDF5.py: globaldb/hdf5 creator
├── losoto
│   ├── operations.py: library of steps losoto can perform
│   └── _version.py
├── losoto.parset: example parset
├── losoto.py: program to run
└── README.md: this file
