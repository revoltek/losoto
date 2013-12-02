LoSoTo: LOFAR solutions tool

Author:
Francesco de Gasperin (overall project)
David Rafferty (overall project)
Reinout Van Weeren (SMOOTH, CLOCKTEC)
Maaijke Mevius (CLOCKTEC)
Bas van der Tol (CLOCKTEC, TOOLS)

This library have tools to deal with LOFAR solution tables.

- set to zero solutions (step: RESET)
- clock/tec separation (step: CLOCKTEC) - TBI
- flag outliers (step: FLAG) - TBI
- solutions smoothing (step: SMOOTH)
- solutions interpolation/rescaling (step: INTERP)
- make some inspection plots (step: PLOT)

Included tools:

- H5parm_creator.py:
    * creates an h5parm file from an MS or a globaldb created with parmdb_collector.py
- H5parm_merge.py:
    * copy a solset from a H5parm files into another one
- parmdb_collector.py:
    * fetches parmdb tables from the cluster

To run it on CEP1:
* A copy of the code is in ~fdg/scripts/losoto/
* Source the tool/lofarinit.csh file which is shipped with the code.
* Prepare a paser starting from the losoto.parset (see also the test.parset)
* check the helps of the scripts in tools/ and of losoto.py
* in case of problems write to fdg@hs.uni-hamburg.de

Used packages:
PyTables version:  3.0.0
HDF5 version:      1.8.4
NumPy version:     1.7.1
Numexpr version:   2.1 (not using Intel's VML/MKL)
Zlib version:      1.2.3.3 (in Python interpreter)
LZO version:       2.06 (Aug 12 2011)
BZIP2 version:     1.0.5 (10-Dec-2007)
Blosc version:     1.2.3 (2013-05-17)
Cython version:    0.19.1
Python version:    2.6.5 (r265:79063, Apr 16 2010, 13:57:41) 
[GCC 4.4.3]
Platform:          linux2-x86_64
Byte-ordering:     little
Detected cores:    8
Default encoding:  ascii
