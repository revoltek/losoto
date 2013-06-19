# Set of operations that LoSoTo can perform
# Each operation is a function

def reset( step, parset, H ):
   """Set specific solutions to 0 (phase, rotang) or 1 (amp)
   """
   npixels = parset.getInt('.'.join(["ExpIon.Steps", step, "npixels"]), 100 )
   extent = parset.getFloat('.'.join(["ExpIon.Steps", step, "extent"]), 0 )


def clocktec( step, parset, H ):
   """Perform clock/tec separation
   """
   raise Exception('Not yet implemented.')


def flag( step, parset, H ):
   """Flag outliers
   """
   raise Exception('Not yet implemented.')


def smooth( step, parset, H ):
   """Smooth solutions
   """
   raise Exception('Not yet implemented.')


def interp( step, parset, H ):
   """Interpolate solutions in freq/time
   """
   raise Exception('Not yet implemented.')


#def write( step, parset, H ):
#   """Write solutions back to an hdf5 file
#   """
#   raise Exception('Not yet implemented.')


def plot( step, parset, H ):
   """Make some inspection plots
   """
   raise Exception('Not yet implemented.')


def apply( step, parset, H ):
   """Apply solution to a specific MS
   """
   raise Exception('Not yet implemented.')
   
   MS = parset.getString('.'.join(["LoSoTo.Steps", step, "MS"]), '' )
   if MS == '':
       print "ERROR: no MS given, cannot correct :("
       sys.exit(1)

   t=table(MSName+'/SPECTRAL_WINDOW/',ack=False)
   chan_freq=t.getcol('CHAN_FREQ')
   Freq_Mean=np.mean(chan_freq)
   t.close()

   print "... reading column %s"%InCol
   t=table(MSName,ack=False)
   data=t.getcol(InCol)
   MSTime=t.getcol("TIME")
   A0=t.getcol("ANTENNA1")
   A1=t.getcol("ANTENNA2")
   t.close()
 
  
