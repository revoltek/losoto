

def operation_settozero( step, parset ):
   """
   """
   npixels = parset.getInt('.'.join(["ExpIon.Steps", step, "npixels"]), 100 )
   extent = parset.getFloat('.'.join(["ExpIon.Steps", step, "extent"]), 0 )


