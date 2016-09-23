LoSoTo: LOFAR solutions tool
============================

Authors:
* Francesco de Gasperin (main developer)
* David Rafferty (h5parm exporter)
* Maaijke Mevius (CT separation)
* Jose Sabater Montes
* Martin Hardcastle
* Andreas Horneffer

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
* Prepare a parset starting from the parset/losoto.parset
* in case of problems write to Francesco de Gasperin: astro@voo.it

Flag:
LoSoTo has an internl flag system encoded in the weight column. Weights are not currently used but some tasks set them to 0 to "flag" a datapoint. Flagged datapoints are then saved as NaNs by the H5parm_exporter.py

Documentation:
Documentation extracted from the LOFAR cookbook is present in the doc directory
A detailed explanation of all the parameters for each operation is in: __parsets/losoto.parset__

Used packages:
* PyTables version:  3.2.2
* HDF5 version:      1.8.4
* NumPy version:     1.9.0
* Python version:    2.6.5 (r265:79063, Apr 16 2010, 13:57:41) [GCC 4.4.3]
* Platform:          linux2-x86_64
* Byte-ordering:     little
* Default encoding:  ascii
