LoSoTo: LOFAR solutions tool
============================

Main authors:
* Francesco de Gasperin (main developer)
* David Rafferty
* Maaijke Mevius

Contents:
* __doc/__: documentation
* __examples/__: some examples h5parm to use with validation/test parsets
* __operations/__: containts all the modules for operations
* __parsets/__: some examples parsets
* __tools/__: contains some tools, mostly to convert back and forward from parmdb to h5parm

Install:
* Get LoSoTo from https://github.com/revoltek/losoto 
* Install losoto by simply: python setup.py install --prefix=~/mydir/
* Alternatively: pip install --allow-external --upgrade --user https://github.com/revoltek/losoto/archive/master.zip 
* In cep3 use the copy of the code in ~fdg/scripts/losoto/ (source the tool/lofarinit.[c]sh file which is shipped with the code)
* Prepare a parset starting from the parset/losoto2.parset
* in case of problems write to Francesco de Gasperin: astro@voo.it

Cite:
* If you use LoSoTo for your data reduction, please acknowledge it with "This work had made use of the Lofar Solution Tool (LoSoTo), developed by F. de Gasperin."

Documentation:
* Documentation of LoSoTo API is at: __http://revoltek.github.io/losoto/losoto.html__
* A detailed explanation of all the parameters is at: __http://revoltek.github.io/losoto/losoto.operations.html__
* An example parset is in: __parsets/losoto2.parset__

Required packages:
* PyTables version:  >3.0.0
* HDF5 version:      >1.8.4
* NumPy version:     >1.9.0
