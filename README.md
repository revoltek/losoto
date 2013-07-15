LoSoTo: LOFAR solutions tool

Author:
Francesco de Gasperin (overall project)
Cyril Tasse (APPLY)
Reinout Van Weeren (SMOOTH)
Maaijke Mevius (CLOCKTEC)
Bas van der Tol (CLOCKTEC)

This library have tools to deal with LOFAR solution tables.

- set to zero solutions (step: RESET)
- clock/tec separation (step: CLOCKTEC)
- flag outliers (step: FLAG)
- solutions smoothing (step: SMOOTH)
- solutions interpolation (step: INTERP)
- make some inspection plots (step: PLOT)
- apply solutions (step: APPLY)

Included tools:

- H5parm_creator.py:
    * fetches parmdb tables from the cluster
    * creates an h5parm file from a parmdb
