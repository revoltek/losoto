LoSoTo: LOFAR solutions tool
============================

### Main developer:
* Francesco de Gasperin

### Contributed code:
* David Rafferty
* Henrik Edler
* Maaijke Mevius
* Jose Sabater Montes
* Martin Hardcastle
* Andreas Horneffer
* and the LOFAR Survey Key Science Project

### Cite:
* If you use LoSoTo for your data reduction, please cite: de Gasperin et al. 2019, Astronomy & Astrophysics, Volume 622, id.A5
(https://ui.adsabs.harvard.edu/abs/2019A&A...622A...5D/abstract)

### Install:
* Install latest official release of LoSoTo from PyPI:
  `pip install --upgrade losoto`
* Or, install latest development version of LoSoTo from GitHub:
  `pip install --upgrade --user https://github.com/revoltek/losoto/archive/master.zip`
* Or clone the LoSoTo repository from https://github.com/revoltek/losoto, and install LoSoTo using:
  `python setup.py install --prefix=~/mydir/`
* Prepare a `parset` file, which can be based on the file `parset/losoto2.parset`
* In case of problems write to Francesco de Gasperin: astro@voo.it

### H5parm plotter (GUI):
This is an external project maintained by Frits Sweijen: https://github.com/tikk3r/lofar-h5plot

### Documentation:
* Documentation of LoSoTo API is at: __http://revoltek.github.io/losoto/losoto.html__
* A detailed explanation of all the parameters is at: __http://revoltek.github.io/losoto/losoto.operations.html__
* A few important articles/how-to are on the github wiki: __https://github.com/revoltek/losoto/wiki__

### Contents:
* __bin/__: the losoto executable and some manipulators for H5parms
* __docs/__: documentation
* __examples/__: some examples h5parm to use with validation/test parsets
* __losoto/operations/__: contains all the modules for operations
* __parsets/__: some examples parsets
* __tools/__: contains some external tools

### Required packages:
* Python 3
* PyTables version:  >3.4.0
* HDF5 version:      >1.8.4
* NumPy version:     >1.9.0
* Scipy              >1.4 (for interpolatedirections)
