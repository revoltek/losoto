# Set of operations that LoSoTo can perform
# Each operation is a function

def reset( step, parset ):
   """Set specific solutions to 0 (phase, rotang) or 1 (amp)
   """
   npixels = parset.getInt('.'.join(["ExpIon.Steps", step, "npixels"]), 100 )
   extent = parset.getFloat('.'.join(["ExpIon.Steps", step, "extent"]), 0 )


def clocktec( step, parset ):
   """Perform clock/tec separation
   """
   raise Exception('Not yet implemented.')

def flag( step, parset ):
   """Flag outliers
   """

def smooth( step, parset ):
   """Smooth solutions
   """

def interp( step, parset ):
   """Interpolate solutions in freq/time
   """

def write( step, parset ):
   """Write solutions back to an hdf5 file
   """

def plot( step, parset ):
   """Make some inspection plots
   """

def apply( step, parset ):
   """Apply solution to a specific MS
   """

