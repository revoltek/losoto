#!/usr/bin/env python
# coding: utf-8

from losoto.h5parm import h5parm
import argparse
import tables

def soltab_swap_freq_time(soltab):
    """Swap the frequency and time axes to make the frequency the fastest varying axis

    Parameters
    ----------
    soltab : Soltab
        Soltab object which will be changed
    """
    vals = soltab.getValues(retAxesVals = False)
    weights = soltab.getValues(weight=True, retAxesVals=False)

    if vals.shape != weights.shape:
        raise RuntimeError("Shape of weights differs from shape of values")

    axesnames = soltab.getAxesNames()
    axesnums = list(range(len(axesnames)))

    if 'freq' not in axesnames or 'time' not in axesnames:
        print("Nothing to be done, no freq + time axes in " + soltab.name)
        return

    freqindex = axesnames.index('freq')
    timeindex = axesnames.index('time')

    if freqindex > timeindex:
        print("Nothing to be done, freq already varies fastest in " + soltab.name)
        return

    # Swap the time and frequency axis in the axes names and numbers
    axesnums[freqindex], axesnums[timeindex] = axesnums[timeindex], axesnums[freqindex]
    axesnames[freqindex], axesnames[timeindex] = axesnames[timeindex], axesnames[freqindex]

    # Swap the axes order in the metadata
    soltab.obj.val._f_setattr("AXES", ",".join(axesnames))
    # Transpose the values
    vals = vals.transpose(axesnums)
    weights = weights.transpose(axesnums)

    # Need to remove the array from the file because changing shape is not supported by pytables
    # Store the attributes in a dict
    attrs = soltab.obj.val._v_attrs
    attrsdict = {}
    for attrname in attrs._f_list():
        attrsdict[attrname] = attrs[attrname]
    soltab.obj.val._f_remove()
    soltab.obj.weight._f_remove()
    # Create new val here
    soltab.obj._v_file.create_array(soltab.obj._v_pathname, 'val', obj=vals, atom=tables.Float64Atom())
    soltab.obj._v_file.create_array(soltab.obj._v_pathname, 'weight', obj=weights, atom=tables.Float16Atom())
    # Restore the original attributes
    for attrname in attrsdict:
        soltab.obj.val._f_setattr(attrname, attrsdict[attrname])

    soltab.addHistory("Swap frequency and time axes to make frequency vary fastest")

    soltab.obj._f_flush()

    print("Successfully swapped frequency and time axes in " + soltab.name)

def h5parm_swap_freq_time(h5parmname, solset='sol000', soltab='all'):
    """Open an H5Parm and swap the frequency and time axes to make the frequency the fastest varying axis

    Parameters
    ----------
    h5parmname : Soltab
        H5Parm filename
    solset : str
        Solset name (default 'sol000')
    soltab : str
        Soltab name or 'all' for all soltabs
    """
    
    h5 = h5parm(h5parmname, False)
    solset = h5.getSolset('sol000')

    if soltab=='all':
        soltabnames = solset.getSoltabNames()
    else:
        soltabnames = [soltab]

    for soltabname in soltabnames:
        soltab = solset.getSoltab(soltabname)
        soltab_swap_freq_time(soltab)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Swap the time and frequency axes in a H5Parm soltab")
    parser.add_argument("h5parmname", help="H5Parm filename")
    parser.add_argument("--solset", "-s", help="Solset to be changed (default 'sol000')", default="sol000")
    parser.add_argument("--soltab", "-t", help="Soltab to be changed, or 'all' for all (default)", default="all")

    args = parser.parse_args()
    h5parm_swap_freq_time(args.h5parmname, solset=args.solset, soltab=args.soltab)
