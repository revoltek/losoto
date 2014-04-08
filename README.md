LoSoTo: LOFAR solutions tool
============================

Authors:
* Francesco de Gasperin
* David Rafferty
* Marting Hardcastle
* Maaijke Mevius
* Bas van der Tol

Contents:
* __doc/__: documentation
* __examples/__: some examples h5parm to use with validation/test parsets
* __operations/__: containts all the modules for operations
* __parsets/__: some examples parsets
* __tools/__: contains some tools, mostly to convert back and forward from parmdb to h5parm

To run it on CEP1:
* A copy of the code is in ~fdg/scripts/losoto/
* Source the tool/lofarinit.csh file which is shipped with the code.
* Prepare a paset starting from the parset/losoto.parset
* check the helps of the scripts in tools/ and of losoto.py
* in case of problems write to fdg@hs.uni-hamburg.de

Used packages:
* PyTables version:  3.0.0
* HDF5 version:      1.8.4
* NumPy version:     1.7.1
* Numexpr version:   2.1 (not using Intel's VML/MKL)
* Zlib version:      1.2.3.3 (in Python interpreter)
* LZO version:       2.06 (Aug 12 2011)
* BZIP2 version:     1.0.5 (10-Dec-2007)
* Blosc version:     1.2.3 (2013-05-17)
* Cython version:    0.19.1
* Python version:    2.6.5 (r265:79063, Apr 16 2010, 13:57:41) [GCC 4.4.3]
* Platform:          linux2-x86_64
* Byte-ordering:     little
* Detected cores:    8
* Default encoding:  ascii
