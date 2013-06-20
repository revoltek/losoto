LoSoTo: LOFAR solutions tool

Author:
Francesco de Gasperin (overall project)
Cyril Tasse (APPLY)
Reinout Van Weeren (SMOOTH)
Maaijke Mevius (CLOCKTEC)
Bas van der Tol (CLOCKTEC)

This library have tools to deal with LOFAR solution tables.

- fetches parmdb tables from the cluster and creates an hdf5 file
- easy access to hdf5 file
- set to zero solutions (step: RESET)
- clock/tec separation (step: CLOCKTEC)
- flag outliers (step: FLAG)
- solutions smoothing (step: SMOOTH)
- solutions interpolation (step: INTERP)
- make some inspection plots (step: PLOT)
- apply solutions (step: APPLY)

.
├── ClassMakeHDF5.py: globaldb/hdf5 creator
├── losoto
│   ├── operations.py: library of steps losoto can perform
│   └── _version.py
├── losoto.parset: example parset
├── losoto.py: program to run
└── README.md: this file
