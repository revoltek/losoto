#!/usr/bin/env python
# coding: utf-8

from losoto.h5parm import h5parm
import tables
import unittest
import numpy as np
import os
from swap_freq_time_axes import h5parm_swap_freq_time

class TestH5parmSwapFreqTime(unittest.TestCase):
    def test_swap_freq_time_axes(self):
      import os, tempfile
      h5fname = tempfile.mktemp(suffix='.h5')

      h5 = h5parm(h5fname, readonly=False)
      solset = h5.makeSolset("sol000")

      dirvals = ["CasA", "VirA"]
      timevals = np.arange(0, 5)
      freqvals = np.arange(0, 3)

      vals = np.ones((2, 3, 5))
      soltab = solset.makeSoltab(soltype=b"phase", soltabName="phase000",
                                 axesNames=["dir","freq","time"], 
                                 axesVals=[dirvals, freqvals, timevals],
                                 vals=vals, weights=vals)

      h5.close()

      h5parm_swap_freq_time(h5fname)

      h5 = h5parm(h5fname, readonly=False)
      solset = h5.getSolset("sol000")
      soltab = solset.getSoltab("phase000")
      vals = soltab.getValues(retAxesVals = False)
      weights = soltab.getValues(retAxesVals=False, weight=True)
      axesnames = soltab.getAxesNames()

      self.assertEqual(vals.shape, (2, 5, 3))
      self.assertEqual(weights.shape, (2, 5, 3))
      self.assertEqual(axesnames, ["dir", "time", "freq"])

      os.remove(h5fname)

if __name__ == '__main__':
    unittest.main()

