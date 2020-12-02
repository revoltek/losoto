#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Retrieving and writing data in H5parm format

import os, sys, re, itertools
import numpy as np
import tables
import losoto._version
from losoto._logging import logger as logging
from losoto.lib_losoto import deprecated_alias
# check for tables version
if int(tables.__version__.split('.')[0]) < 3:
    logging.critical('pyTables version must be >= 3.0.0, found: '+tables.__version__)
    sys.exit(1)


def openSoltab(h5parmFile, solsetName=None, soltabName=None, address=None, readonly=True):
    """
    Convenience function to get a soltab object from an h5parm file and an address like "solset000/phase000".

    Parameters
    ----------
    h5parmFile : str
        H5parm filename.
    solsetName : str
        solset name
    soltabName : str
        soltab name
    address : str
        solset/soltab name (to use in place of the parameters solset and soltab).
    readonly : bool, optional
        if True the table is open in readonly mode, by default True.

    Returns
    -------
    Soltab obj
        A solution table object.
    """
    h5 = h5parm(h5parmFile, readonly)
    if solsetName is None or soltabName is None:
        if address is None:
            logging.error('Address must be specified if solsetName and soltabName are not given.')
            sys.exit(1)
        solsetName, soltabName = address.split('/')
    solset = h5.getSolset(solsetName)
    return solset.getSoltab(soltabName)


class h5parm( object ):
    """
    Create an h5parm object.

    Parameters
    ----------
    h5parmFile : str
        H5parm filename.
    readonly : bool, optional
        if True the table is open in readonly mode, by default True.
    complevel : int, optional
        compression level from 0 to 9 when creating the file, by default 5.
    complib : str, optional
        library for compression: lzo, zlib, bzip2, by default zlib.
    """

    def __init__(self, h5parmFile, readonly=True, complevel=0, complib='zlib'):

        self.H = None  # variable to store the pytable object
        self.fileName = h5parmFile

        if os.path.isfile(h5parmFile):
            if not tables.is_hdf5_file(h5parmFile):
                logging.critical('Not a HDF5 file: '+h5parmFile+'.')
                raise Exception('Not a HDF5 file: '+h5parmFile+'.')
            if readonly:
                logging.debug('Reading from '+h5parmFile+'.')
                self.H = tables.open_file(h5parmFile, 'r', IO_BUFFER_SIZE=1024*1024*10, BUFFER_TIMES=500)
            else:
                logging.debug('Appending to '+h5parmFile+'.')
                self.H = tables.open_file(h5parmFile, 'r+', IO_BUFFER_SIZE=1024*1024*10, BUFFER_TIMES=500)

            # Check if it's a valid H5parm file: attribute h5parm_version should be defined in any solset
            is_h5parm = True
            for node in self.H.root:
                if 'h5parm_version' not in node._v_attrs:
                    is_h5parm=False
                    break
            if not is_h5parm:
                logging.warning('Missing H5pram version. Is this a properly made h5parm?')

        else:
            if readonly:
                raise Exception('Missing file '+h5parmFile+'.')
            else:
                logging.debug('Creating '+h5parmFile+'.')
                # add a compression filter
                f = tables.Filters(complevel=complevel, complib=complib)
                self.H = tables.open_file(h5parmFile, filters=f, mode='w', IO_BUFFER_SIZE=1024*1024*10, BUFFER_TIMES=500)


    def close(self):
        """
        Close the open table.
        """
        logging.debug('Closing table.')
        self.H.close()


    def __str__(self):
        """
        Returns
        -------
        string
            Info about H5parm contents.
        """
        return self.printInfo()


    def makeSolset(self, solsetName=None, addTables=True):
        """
        Create a new solset, if the provided name is not given or exists
        then it falls back on the first available sol###.

        Parameters
        ----------
        solset : str
            Name of the solution set.
        addTables : bool, optional
            If True add antenna/direction/array tables, by default True.

        Returns
        -------
        solset obj
            Newly created solset object.
        """

        if type(solsetName) is str and not re.match(r'^[A-Za-z0-9_-]+$', solsetName):
            logging.warning('Solution-set '+solsetName+' contains unsuported characters. Use [A-Za-z0-9_-]. Switching to default.')
            solsetName = None

        if solsetName in self.getSolsetNames():
            logging.warning('Solution-set '+solsetName+' already present. Switching to default.')
            solsetName = None

        if solsetName is None:
            solsetName = self._firstAvailSolsetName()

        logging.info('Creating a new solution-set: '+solsetName+'.')
        solset = self.H.create_group("/", solsetName)
        solset._f_setattr('h5parm_version', losoto._version.__h5parmVersion__)

        if addTables:
            # add antenna table
            logging.info('--Creating new antenna table.')
            descriptor = np.dtype([('name', np.str_, 16),('position', np.float32, 3)])
            soltab = self.H.create_table(solset, 'antenna', descriptor, \
                    title = 'Antenna names and positions', expectedrows = 50)

            # add direction table
            logging.info('--Creating new source table.')
            descriptor = np.dtype([('name', np.str_, 128),('dir', np.float32, 2)])
            soltab = self.H.create_table(solset, 'source', descriptor, \
                    title = 'Source names and directions', expectedrows = 25)

        return Solset(solset)


    def getSolsets(self):
        """
        Get all solution set objects.

        Returns
        -------
        list
            A list of all solsets objects.
        """
        solsets = []
        for solset in self.H.root._v_groups.values():
            solsets.append(Solset(solset))
        return solsets


    def getSolsetNames(self):
        """
        Get all solution set names.

        Returns
        -------
        list
            A list of str of all solsets names.
        """
        solsetNames = []
        for solsetName in iter( list(self.H.root._v_groups.keys()) ):
            solsetNames.append(solsetName)
        return solsetNames


    def getSolset(self, solset):
        """
        Get a solution set with a given name.

        Parameters
        ----------
        solset : str
            Name of the solution set.

        Returns
        -------
        solset obj
            Return solset object.
        """
        if not solset in self.getSolsetNames():
            logging.critical("Cannot find solset: "+solset+".")
            raise Exception("Cannot find solset: "+solset+".")

        return Solset(self.H.get_node('/',solset))


    def _firstAvailSolsetName(self):
        """
        Find the first available solset name which has the form of "sol###".

        Returns
        -------
        str
            Solset name.
        """
        nums = []
        for solsetName in self.getSolsetNames():
            if re.match(r'^sol[0-9][0-9][0-9]$', solsetName):
                nums.append(int(solsetName[-3:]))

        return "sol%03d" % min(list(set(range(1000)) - set(nums)))


    def printInfo(self, filter=None, verbose=False):
        """
        Used to get readable information on the h5parm file.

        Parameters
        ----------
        filter: str, optional
            Solution set name to get info for
        verbose: bool, optional
            If True, return additional info on axes

        Returns
        -------
        str
            Returns a string with info about H5parm contents.
        """
        if (sys.version_info > (3, 0)):
            from itertools import zip_longest
        else:
            from itertools import izip_longest as zip_longest

        def grouper(n, iterable, fillvalue=' '):
            """
            Groups iterables into specified groups

            Parameters
            ----------
            n : int
                number of iterables to group
            iterable : iterable
                iterable to group
            fillvalue : str
                value to use when to fill blanks in output groups

            Example
            -------
            grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
            """
            args = [iter(iterable)] * n
            return zip_longest(fillvalue=fillvalue, *args)

        def wrap(text, width=80):
            """
            Wraps text to given width and returns list of lines
            """
            lines = []
            for paragraph in text.split('\n'):
                line = []
                len_line = 0
                for word in paragraph.split(' '):
                    word.strip()
                    len_word = len(word)
                    if len_line + len_word <= width:
                        line.append(word)
                        len_line += len_word + 1
                    else:
                        lines.append(' '.join(line))
                        line = [21*' '+word]
                        len_line = len_word + 22
                lines.append(' '.join(line))
            return lines

        info = "\nSummary of %s\n" % self.fileName
        solsets = self.getSolsets()

        # Filter on solset name
        if filter is not None:
            keys_to_remove = []
            info += "\nFiltering on solution set name with filter = '{0}'\n".format(filter)
            solsets = [solset for solset in solsets if re.search(filter, solset.name)]

        if len(solsets) == 0:
            info += "\nNo solution sets found.\n"
            return info

        # delete axes value file if already present
        if verbose and os.path.exists(self.fileName+'-axes_values.txt'):
            logging.warning('Overwriting '+self.fileName+'-axes_values.txt')
            os.system('rm '+self.fileName+'-axes_values.txt')

        # For each solution set, list solution tables, sources, and antennas
        for solset in solsets:
            info += "\nSolution set '%s':\n" % solset.name
            info += "=" * len(solset.name) + "=" * 16 + "\n\n"

            # Add direction (source) names
            sources = sorted( solset.getSou().keys() )
            info += "Directions: "
            for src_name1, src_name2, src_name3 in grouper(3, sources):
                info += "{0:}\t{1:}\t{2:}\n            ".format(src_name1, src_name2, src_name3)

            # Add station names
            antennas = sorted( solset.getAnt().keys() )
            info += "\nStations: "
            for ant1, ant2, ant3, ant4 in grouper(4, antennas):
                info += "{0:}\t{1:}\t{2:}\t{3:}\n          ".format(ant1, ant2, ant3, ant4)

            # For each table, add length of each axis and history of
            # operations applied to the table.
            if verbose:
                logging.warning('Axes values saved in '+self.fileName+'-axes_values.txt')
                f = open(self.fileName+'-axes_values.txt','a')
            soltabs = solset.getSoltabs()
            names = np.array([s.name for s in soltabs])
            sorted_soltabs = [x for _, x in sorted(zip(names, soltabs))]
            for soltab in sorted_soltabs:
                try:
                    if verbose:
                        f.write("### /"+solset.name+"/"+soltab.name+"\n")
                    logging.debug('Fetching info for '+soltab.name+'.')
                    axisNames = soltab.getAxesNames()
                    axis_str_list = []
                    for axisName in axisNames:
                        nslots = soltab.getAxisLen(axisName)
                        if nslots > 1:
                            pls = "s"
                        else:
                            pls = ""
                        axis_str_list.append("%i %s%s" % (nslots, axisName, pls))
                        if verbose:
                            f.write(axisName+": ")
                            vals = soltab.getAxisValues(axisName)
                            # ugly hardcoded workaround to print all the important decimal values for time/freq
                            if axisName == 'freq': f.write(" ".join(["{0:.8f}".format(v) for v in vals])+"\n\n")
                            elif axisName == 'time': f.write(" ".join(["{0:.7f}".format(v) for v in vals])+"\n\n")
                            else: f.write(" ".join(["{}".format(v) for v in vals])+"\n\n")
                    info += "\nSolution table '%s' (type: %s): %s\n" % (soltab.name, soltab.getType(), ", ".join(axis_str_list))
                    weights = soltab.getValues(weight = True, retAxesVals = False)
                    vals = soltab.getValues(weight = False, retAxesVals = False)
                    info += '    Flagged data: %.3f%%\n' % (100.*np.sum(weights==0 | np.isnan(vals))/len(weights.flat))

                    # Add some extra attributes stored in screen-type tables
                    if soltab.getType() == 'screen':
                        attr_names = soltab.obj._v_attrs._v_attrnames
                        add_head = True
                        for n in attr_names:
                            if n in ['beta', 'freq', 'height', 'order']:
                                if add_head:
                                    info += '    Screen attributes:\n'
                                    add_head = False
                                info += '        {0}: {1}\n'.format(n, soltab.obj._v_attrs[n])

                    # Add history
                    history = soltab.getHistory()
                    if history != "":
                        info += 4*" " + "History: "
                        joinstr = "\n" + 13*" "
                        info += joinstr.join(wrap(history)) + "\n"
                except tables.exceptions.NoSuchNodeError:
                    info += "\nSolution table '%s': No valid data found\n" % (soltab.name)

            if verbose:
                f.close()
        return info


class Solset( object ):
    """
    Create a solset object

    Parameters
    ----------
    solset : pytables group
        The solution set pytables group object
    """

    def __init__(self, solset):

        if not isinstance( solset, tables.Group ):
            logging.error("Object must be initialized with a pyTables Group object.")
            sys.exit(1)

        self.obj = solset
        self.name = solset._v_name

    def close(self):
        """
        """
        self.obj._g_flushGroup()


    def delete(self):
        """
        Delete this solset.
        """
        logging.info("Solset \""+self.name+"\" deleted.")
        self.obj._f_remove(recursive=True)


    def rename(self, newname, overwrite=False):
        """
        Rename this solset.

        Parameters
        ----------
        newname : str
            New solution set name.
        overwrite : bool, optional
            Overwrite existing solset with same name.
        """
        self.obj._f_rename(newname, overwrite)
        logging.info('Solset "'+self.name+'" renamed to "'+newname+'".')
        self.name = self.obj._v_name


    def makeSoltab(self, soltype=None, soltabName=None,
            axesNames = [], axesVals = [], chunkShape=None, vals=None,
            weights=None, parmdbType='', weightDtype='f16'):
        """
        Create a Soltab into this solset.

        Parameters
        ----------
        soltype : str
            Solution-type (e.g. amplitude, phase)
        soltabName : str, optional
            The solution-table name, if not specified is generated from the solution-type
        axesNames : list
            List with the axes names
        axesVals : list
            List with the axes values (each is a separate np.array)
        chunkShape : list, optional
            List with the chunk shape
        vals : numpy array
            Array with shape given by the axesVals lenghts
        weights : numpy array
            Same shape of the vals array
            0->FLAGGED, 1->MAX_WEIGHT
        parmdbType : str
            Original parmdb solution type
        weightDtype : str
            THe dtype of weights allowed values are ('f16' or 'f32' or 'f64')

        Returns
        -------
        soltab obj
            Newly created soltab object
        """

        if soltype is None:
            raise Exception("Solution-type not specified while adding a solution-table.")

        # checks on the soltab
        if type(soltabName) is str and not re.match(r'^[A-Za-z0-9_-]+$', soltabName):
            logging.warning('Solution-table '+soltabName+' contains unsuported characters. Use [A-Za-z0-9_-]. Switching to default.')
            soltabName = None

        if soltabName in self.getSoltabNames():
            logging.warning('Solution-table '+soltabName+' already present. Switching to default.')
            soltabName = None

        if soltabName is None:
            soltabName = self._fisrtAvailSoltabName(soltype)

        logging.info('Creating a new solution-table: '+soltabName+'.')

        # check input
        assert len(axesNames) == len(axesVals)
        dim = []
        for i, axisName in enumerate(axesNames):
            dim.append(len(axesVals[i]))
        assert dim == list(vals.shape)
        assert dim == list(weights.shape)

        # if input is OK, create table
        soltab = self.obj._v_file.create_group("/"+self.name, soltabName, title=soltype)
        soltab._v_attrs['parmdb_type'] = parmdbType
        for i, axisName in enumerate(axesNames):
            axis = self.obj._v_file.create_array('/'+self.name+'/'+soltabName, axisName, obj=np.array(axesVals[i]))

        # create the val/weight Carrays
        #val = self.obj._v_file.create_carray('/'+self.name+'/'+soltabName, 'val', obj=vals.astype(np.float64), chunkshape=None, atom=tables.Float64Atom())
        #weight = self.obj._v_file.create_carray('/'+self.name+'/'+soltabName, 'weight', obj=weights.astype(np.float16), chunkshape=None, atom=tables.Float16Atom())
        # array do not have compression but are much faster
        val = self.obj._v_file.create_array('/'+self.name+'/'+soltabName, 'val', obj=vals.astype(np.float64), atom=tables.Float64Atom())
        assert weightDtype in ['f16','f32', 'f64'], "Allowed weight dtypes are 'f16','f32', 'f64'"
        if weightDtype == 'f16':
            np_d = np.float16
            pt_d = tables.Float16Atom()
        elif weightDtype == 'f32':
            np_d = np.float32
            pt_d = tables.Float32Atom()
        elif weightDtype == 'f64':
            np_d = np.float64
            pt_d = tables.Float64Atom()
        weight = self.obj._v_file.create_array('/'+self.name+'/'+soltabName, 'weight', obj=weights.astype(np_d), atom=pt_d)
        axisNames = ','.join([axisName for axisName in axesNames])
        val.attrs['AXES'] = axisNames.encode()
        weight.attrs['AXES'] = axisNames.encode()

        return Soltab(soltab)


    def _fisrtAvailSoltabName(self, soltype):
        """
        Find the first available soltab name which
        has the form of "soltypeName###"

        Parameters
        ----------
        soltype : str
            Type of solution (amplitude, phase, RM, clock...)

        Returns
        -------
        str
            First available soltab name
        """
        nums = []
        for soltab in self.getSoltabs():
            if re.match(r'^'+soltype+'[0-9][0-9][0-9]$', soltab.name):
                nums.append(int(soltab.name[-3:]))

        return soltype+"%03d" % min(list(set(range(1000)) - set(nums)))


    def getSoltabs(self, useCache=False, sel={}):
        """
        Get all Soltabs in this Solset.

        Parameters
        ----------
        useCache : bool, optional
            soltabs obj will use cache, by default False
        sel : dict, optional
            selection dict, by default no selection

        Returns
        -------
        list
            List of solution tables objects for all available soltabs in this solset
        """
        soltabs = []
        for soltab in self.obj._v_groups.values():
            soltabs.append(Soltab(soltab, useCache, sel))
        return soltabs


    def getSoltabNames(self):
        """
        Get all Soltab names in this Solset.

        Returns
        -------
        list
            List of str for all available soltabs in this solset.
        """
        soltabNames = []
        for soltabName in iter( list(self.obj._v_groups.keys()) ):
            soltabNames.append(soltabName)
        return soltabNames


    def getSoltab(self, soltab, useCache=False, sel={}):
        """
        Get a soltab with a given name.

        Parameters
        ----------
        soltab : str
            A solution table name.
        useCache : bool, optional
            Soltabs obj will use cache, by default False.
        sel : dict, optional
            Selection dict, by default no selection.

        Returns
        -------
        soltab obj
            A solution table obj.
        """
        if soltab is None:
            raise Exception("Solution-table not specified while querying for solution-table.")

        if not soltab in self.getSoltabNames():
            raise Exception("Solution-table "+soltab+" not found in solset "+self.name+".")

        return Soltab(self.obj._f_get_child(soltab), useCache, sel)


    def getAnt(self):
        """
        Get the antenna subtable with antenna names and positions.

        Returns
        -------
        dict
            Available antennas in the form {name1:[position coords], name2:[position coords], ...}.
        """
        ants = {}
        try:
            for x in self.obj.antenna:
                ants[x['name'].decode()] = x['position']
        except: pass

        return ants


    def getSou(self):
        """
        Get the source subtable with direction names and coordinates.

        Returns
        -------
        dict
            Available sources in the form {name1:[ra,dec], name2:[ra,dec], ...}.
        """
        sources = {}
        try:
            for x in self.obj.source:
                sources[x['name'].decode()] = x['dir']
        except: pass

        return sources

    def getAntDist(self, ant=None):
        """
        Get antenna distance to a specified one.

        Parameters
        ----------
        ant : str
            An antenna name.

        Returns
        -------
        str
            Dict of distances to each antenna. The distance with the antenna "ant" is 0.
        """
        if ant is None:
            raise "Missing antenna name."

        ants = self.getAnt()

        if not ant in list(ants.keys()):
            raise "Missing antenna %s in antenna table." % ant 

        return {a:np.sqrt( (loc[0]-ants[ant][0])**2 + (loc[1]-ants[ant][1])**2 + (loc[2]-ants[ant][2])**2 ) for a, loc in ants.items() }



class Soltab( object ):
    """
    Parameters
    ----------
    soltab : pytables Table obj
        Pytable Table object.
    useCache : bool, optional
        Cache all data in memory, by default False.
    **args : optional
        Used to create a selection.
        Selections examples:
        axisName = None # to select all
        axisName = xxx # to select ONLY that value for an axis
        axisName = [xxx, yyy, zzz] # to selct ONLY those values for an axis
        axisName = 'xxx' # regular expression selection
        axisName = {min: xxx} # to selct values grater or equal than xxx
        axisName = {max: yyy} # to selct values lower or equal than yyy
        axisName = {min: xxx, max: yyy} # to selct values greater or equal than xxx and lower or equal than yyy
    """

    def __init__(self, soltab, useCache = False, args = {}):

        if not isinstance( soltab, tables.Group ):
            logging.error("Object must be initialized with a pyTables Table object.")
            sys.exit(1)

        self.obj = soltab
        self.name = soltab._v_name

        # list of axes names, set once to speed up calls
        axesNamesInH5 = soltab.val.attrs['AXES'].decode()
        self.axesNames = axesNamesInH5.split(',')

        # dict of axes values, set once to speed up calls (a bit of memory usage though)
        self.axes = {}
        for axis in self.getAxesNames():
            self.axes[axis] = soltab._f_get_child(axis)

        # initialize selection
        self.setSelection(**args)

        self.useCache = useCache
        if self.useCache:
            logging.debug("Caching...")
            self.setCache(self.obj.val, self.obj.weight)

        self.fullyFlaggedAnts = None # this is populated if required by reference


    def delete(self):
        """
        Delete this soltab.
        """
        logging.info("Soltab \""+self.name+"\" deleted.")
        self.obj._f_remove(recursive=True)


    def rename(self, newname, overwrite=False):
        """
        Rename this soltab.

        Parameters
        ----------
        newname : str
            New solution table name.
        overwrite : bool, optional
            Overwrite existing soltab with same name.
        """
        self.obj._f_rename(newname, overwrite)
        logging.info('Soltab "'+self.name+'" renamed to "'+newname+'".')
        self.name = self.obj._v_name


    def setCache(self, val, weight):
        """
        Set cache values.

        Parameters
        ----------
        val : array
        weight : array
        """
        self.cacheVal = np.copy(val)
        self.cacheWeight = np.copy(weight)


    def getSolset(self):
        """
        This is used to obtain the parent solset object to e.g. get antennas or create new soltabs.

        Returns
        -------
        solset obj
            A solset obj.
        """
        return Solset(self.obj._v_parent)


    def getAddress(self):
        """
        Get the "solset000/soltab000" type string for this Soltab.

        Returns
        -------
        str
            The solset/soltab address of self.obj as a string.
        """
        return self.obj._v_pathname[1:]


    def clearSelection(self):
        """
        Clear selection, all values are now considered.
        """
        self.setSelection()


    def setSelection(self, update=False, **args):
        """
        Set a selection criteria. For each axes there can be a:
            * string: regexp
            * list: use only the listed values
            * dict: with min/max/[step] to select a range.

        Parameters
        ----------
        **args :
            Valid axes names of the form: pol='XX', ant=['CS001HBA','CS002HBA'], time={'min':1234,'max':'2345','step':4}.
        
        update : bool
            Only update axes passed as arguments, the rest is maintained. Default: False.
            
        """
        # create an initial selection which selects all values
        if not update:
            self.selection = [slice(None)] * len(self.getAxesNames())

        for axis, selVal in iter( list(args.items()) ):
            # if None continue and keep all the values
            if selVal is None: continue
            if not axis in self.getAxesNames():
                logging.warning("Cannot select on axis "+axis+", it doesn't exist. Ignored.")
                continue

            # find the index of the working axis
            idx = self.getAxesNames().index(axis)

            # slice -> let the slice be as it is
            if isinstance(selVal, slice):
                self.selection[idx] = selVal

            # string -> regular expression
            elif type(selVal) is str:
                if not self.getAxisType(axis).char == 'S':
                    logging.warning("Cannot select on axis \""+axis+"\" with a regular expression. Use all available values.")
                    continue
                self.selection[idx] = [i for i, item in enumerate(self.getAxisValues(axis)) if re.search(selVal, item)]

                # transform list of 1 element in a relative slice(), faster as it gets reference
                if len(self.selection[idx]) == 1: self.selection[idx] = slice(self.selection[idx][0],self.selection[idx][0]+1)

            # dict -> min max
            elif type(selVal) is dict:
                axisVals = self.getAxisValues(axis)
                # some checks
                if 'min' in selVal and selVal['min'] > np.max(axisVals):
                    logging.error("Selection with min > than maximum value. Use all available values.")
                    continue
                if 'max' in selVal and selVal['max'] < np.min(axisVals):
                    logging.error("Selection with max < than minimum value. Use all available values.")
                    continue

                if 'min' in selVal and 'max' in selVal:
                    self.selection[idx] = slice(np.where(axisVals >= selVal['min'])[0][0], np.where(axisVals <= selVal['max'])[0][-1]+1)
                    #thisSelection[idx] = list(np.where((axisVals>=selVal['min']) & (axisVals<=selVal['max']))[0])
                elif 'min' in selVal:
                    self.selection[idx] = slice(np.where(axisVals >= selVal['min'])[0][0], None)
                    #thisSelection[idx] = list(np.where(axisVals>=selVal['min'])[0])
                elif 'max' in selVal:
                    self.selection[idx] = slice(0, np.where(axisVals <= selVal['max'])[0][-1]+1)
                    #thisSelection[idx] = list(np.where(axisVals<=selVal['max'])[0])
                else:
                    logging.error("Selection with a dict must have 'min' and/or 'max' entry. Use all available values.")
                    continue
                if 'step' in selVal:
                    self.selection[idx] = slice(self.selection[idx].start, self.selection[idx].stop, selVal['step'])
                    #thisSelection[idx] = thisSelection[idx][::selVal['step']]

            # single val/list -> exact matching
            else:
                if type(selVal) is np.array or type(selVal) is np.ndarray: selVal = selVal.tolist()
                if not type(selVal) is list: selVal = [selVal]
                # convert to correct data type (from parset everything is a str)
                if not self.getAxisType(axis).type is np.string_:
                    selVal = np.array(selVal, dtype=self.getAxisType(axis))
                else:
                    selVal = np.array(selVal)

                if len(selVal) == 1:
                    # speedup in the common case of a single value
                    if not selVal[0] in self.getAxisValues(axis).tolist():
                        logging.error('Cannot find value %s in axis %s. Skip selection.' % (selVal[0], axis))
                        return
                    self.selection[idx] = [self.getAxisValues(axis).tolist().index(selVal)]
                else:
                    self.selection[idx] = [i for i, item in enumerate(self.getAxisValues(axis)) if item in selVal]

                # transform list of 1 element in a relative slice(), faster as it gets a reference
                if len(self.selection[idx]) == 1: self.selection[idx] = slice(self.selection[idx][0], self.selection[idx][0]+1)
                # transform list of continuous numbers in slices, faster as it gets a reference
                elif len(self.selection[idx]) != 0 and len(self.selection[idx])-1 == self.selection[idx][-1] - self.selection[idx][0]:
                    self.selection[idx] = slice(self.selection[idx][0], self.selection[idx][-1]+1)

            # if a selection return an empty list (maybe because of a wrong name), then use all values
            if type(self.selection[idx]) is list and len(self.selection[idx]) == 0:
                logging.warning("Empty/wrong selection on axis \""+axis+"\". Use all available values.")
                self.selection[idx] = slice(None)


    def getType(self):
        """
        Get the solution type of this Soltab.

        Returns
        -------
        str
            Return the type of the solution-tables (e.g. amplitude).
        """

        return self.obj._v_title


    def getAxesNames(self):
        """
        Get axes names.

        Returns
        -------
        list
            A list with all the axis names in the correct order for
            slicing the getValuesGrid() reurned matrix.
        """

        return self.axesNames[:]


    def getAxisLen(self, axis, ignoreSelection=False):
        """
        Return an axis lenght.

        Parameters
        ----------
        axis : str
            The name of the axis.
        ignoreSelection : bool, optional
            If True returns the axis lenght without any selection active, by default False.

        Returns
        -------
        int
            The axis lenght.
        """
        return len(self.getAxisValues(axis, ignoreSelection = ignoreSelection))


    def getAxisType(self, axis):
        """
        Return the axis dtype

        Parameters
        ----------
        axis : str
            The name of the axis.

        Returns
        -------
        dtype
            The axis dtype.
        """
        if axis not in self.getAxesNames():
            logging.error('Axis \"'+axis+'\" not found.')
            return None

        return self.obj._f_get_child(axis).dtype


    def getAxisValues(self, axis, ignoreSelection=False):
        """
        Get the values of a given axis.

        Parameters
        ----------
        axis : str
            The name of the axis.
        ignoreSelection : bool, optional
            If True returns the axis values without any selection active, by default False.

        Returns
        -------
        list
            A copy of all values present along a specific axis.
        """
        if axis not in self.getAxesNames():
            logging.error('Axis \"'+axis+'\" not found.')
            return None

        if ignoreSelection:
            axisvalues = np.copy(self.axes[axis])
        else:
            axisIdx = self.getAxesNames().index(axis)
            axisvalues = np.copy(self.axes[axis][ self.selection[axisIdx] ])

        if axisvalues.dtype.str[0:2] == '|S':
            # Convert to native string format for python 3
            return axisvalues.astype(str)
        else:
            return axisvalues


    def setAxisValues(self, axis, vals):
        """
        Set the value of a specific axis

        Parameters
        ----------
        axis : str
            The name of the axis.
        vals : array
            Values
        """

        if axis not in self.getAxesNames():
            logging.error('Axis \"'+axis+'\" not found.')
            return

        axisIdx = self.getAxesNames().index(axis)
        self.axes[axis][ self.selection[axisIdx] ] = vals


    def setValues(self, vals, selection = None, weight = False):
        """
        Save values in the val grid

        Parameters
        ----------
        vals : array, float
            values to write as an n-dimentional array which match the selection dimention
            if a float is passed or the selected data are set to that value

        selection : selection format, optional
            To set only a subset of data, overriding global selectioan, by default use global selection.
            This is used to set values in a loop of getValueIter(). Global seclection is NOT overwritten.

        weight : bool, optional
            If true store in the weights instead that in the vals, by default False
        """
        if selection is None: selection = self.selection

        if self.useCache:
            if weight: dataVals = self.cacheWeight
            else: dataVals = self.cacheVal
        else:
            if weight: dataVals = self.obj.weight
            else: dataVals = self.obj.val

        # NOTE: pytables has a nasty limitation that only one list can be applied when selecting.
        # Conversely, one can apply how many slices he wants.
        # Single values/contigous values are converted in slices in h5parm.
        # This try/except implements a workaround for this limitation. Once the pytables will be updated, the except can be removed.
        try:
            # the float check allows quick reset of large arrays to a single value
            if isinstance(vals, (np.floating, float)): dataVals[tuple(selection)] = vals
            # the reshape is needed when saving e.g. [512] (vals shape) into [512,1,1] (selection output)
            else: dataVals[tuple(selection)] = np.reshape(vals, dataVals[tuple(selection)].shape)
        except:
            #logging.debug('Optimizing selection writing '+str(selection))
            selectionListsIdx = [i for i, s in enumerate(selection) if type(s) is list]
            subSelection = selection[:]
            # create a subSelection also for the "vals" array
            subSelectionForVals = [slice(None) for i in range(len(subSelection))]
            # cycle across lists and save data index by index
            for selectionListValsIter in itertools.product(*[selection[selectionListIdx] for selectionListIdx in selectionListsIdx[1:]]):
                for i, selectionListIdx in enumerate(selectionListsIdx[1:]):
                    # this is the sub selection which has a slice for every slice and a single value for every list
                    subSelection[selectionListIdx] = selectionListValsIter[i]
                    subSelectionForVals[selectionListIdx] = i
                if type(vals) == float: dataVals[tuple(subSelection)] = vals
                else: dataVals[tuple(subSelection)] = vals[tuple(subSelectionForVals)]

    def flush(self):
        """
        Copy cached values into the table
        """
        if not self.useCache:
            logging.error("Flushing non cached data.")
            sys.exit(1)

        logging.info("Writing results...")
        self.obj.weight[:] = self.cacheWeight
        self.obj.val[:] = self.cacheVal


    def __getattr__(self, axis):
        """
        Links any attribute with an "axis name" to getValuesAxis("axis name")
        also links val and weight to the relative arrays.

        Parameters
        ----------
        axis : str
            The axis name.
        """
        if axis == 'val':
            return self.getValues(retAxesVals=False)
        elif axis == 'weight':
            return self.getValues(retAxesVals=False, weight=True)
        elif axis in self.getAxesNames():
            return self.getAxisValues(axis)
        else:
            return object.__getattribute__(self, axis)
            #logging.error("Cannot find axis \""+axis+"\".")
            #return None


    def _applyAdvSelection(self, data, selection):
        # NOTE: pytables has a nasty limitation that only one list can be applied when selecting.
        # Conversely, one can apply how many slices he wants.
        # Single values/contigous values are converted in slices in h5parm.
        # This implements a workaround for this limitation. Once the pytables will be updated, the except can be removed.
        if np.sum( [1. for sel in selection if type(sel) is list] ) > 1 and \
           ( type(data) is np.ndarray or \
           np.sum( [len(sel)-1 for sel in selection if type(sel) is list] ) > 0 ):
        
            #logging.debug('Optimizing selection reading '+str(selection))
            # for performances is important to minimize the fetched data
            # move all slices at the first selection and lists afterwards (first list is allowd in firstselection)
            selectionListsIdx = [i for i, s in enumerate(selection) if type(s) is list]
            firstSelection = selection[:]
            for i in selectionListsIdx[1:]:
                firstSelection[i] = slice(None)
            # create a second selection using np.ix_
            secondSelection = []
            for i, sel in enumerate(selection):
                #if i == selectionListsIdx[0]: secondSelection.append(range(self.getAxisLen(self.getAxesNames()[i], ignoreSelection=False)))
                if i == selectionListsIdx[0]: secondSelection.append(list(range(len(sel))))
                elif type(sel) is list: secondSelection.append(sel)
                elif type(sel) is slice: secondSelection.append(list(range(self.getAxisLen(self.getAxesNames()[i], ignoreSelection=False))))
            #print firstSelection
            #print secondSelection
            #print data[tuple(firstSelection)].shape
            #print data[tuple(firstSelection)][np.ix_(*secondSelection)].shape
            return data[tuple(firstSelection)][np.ix_(*secondSelection)]
        else:
            return data[tuple(selection)]


    def _getFullyFlaggedAnts(self):
        if self.fullyFlaggedAnts is None:
            self.fullyFlaggedAnts = [] # fully flagged antennas
            antAxis = self.getAxesNames().index('ant')

            if self.useCache:
                dataWeights = self.cacheWeight
            else:
                dataWeights = self.obj.weight

            for antToCheck in self.getAxisValues('ant', ignoreSelection=True):
                # fully flagged?
                refSelection = [slice(None)] * len(self.getAxesNames())
                refSelection[antAxis] = [self.getAxisValues('ant', ignoreSelection=True).tolist().index(antToCheck)]
                if (self._applyAdvSelection(dataWeights, refSelection) == 0 ).all():
                    self.fullyFlaggedAnts.append(antToCheck)

        return self.fullyFlaggedAnts


    @deprecated_alias(reference='refAnt') # Add alias for backwards compatibility
    def getValues(self, retAxesVals=True, weight=False, refAnt=None, refDir=None):
        """
        Creates a simple matrix of values. Fetching a copy of all selected rows into memory.

        Parameters
        ----------
        retAxesVals : bool, optional
            If true returns also the axes vals as a dict of:
            {'axisname1':[axisvals1],'axisname2':[axisvals2],...}.
            By default True.
        weight : bool, optional
            If true get the weights instead that the vals, by defaul False.
        refAnt : str, optional
            In case of phase or rotation solutions, reference to this station name. By default no reference.
            If "closest" reference to the closest antenna.
        refDir : str, optional
            In case of phase or rotation solutions, reference to this Direction. By default no reference.
            If "center", reference to the central direction.

        Returns
        -------
        array
            A numpy ndarray (values or weights depending on parameters)
            If selected, returns also the axes values
        """
        if self.useCache:
            if weight: dataVals = self.cacheWeight
            else: dataVals = self.cacheVal
        else:
            if weight: dataVals = self.obj.weight
            else: dataVals = self.obj.val

        dataVals = self._applyAdvSelection(dataVals, self.selection)

        # CASE 1: Reference only to ant
        if refAnt and not refDir:
            # TODO: Should there be a warning if only ant is referenced but multiple directions are present?
            if not self.getType() in ['phase', 'scalarphase', 'rotation', 'tec', 'clock', 'tec3rd', 'rotationmeasure']:
                logging.error('Reference possible only for phase, scalarphase, clock, tec, tec3rd, rotation and rotationmeasure solution tables. Ignore referencing.')
            elif not 'ant' in self.getAxesNames():
                logging.error('Cannot find antenna axis for referencing phases. Ignore referencing.')
            elif not refAnt in self.getAxisValues('ant', ignoreSelection = True) and refAnt != 'closest':
                logging.error('Cannot find antenna '+refAnt+'. Ignore referencing.')
            else:
                if self.useCache:
                    if weight: dataValsRef = self.cacheWeight
                    else: dataValsRef = self.cacheVal
                else:
                    if weight: dataValsRef = self.obj.weight
                    else: dataValsRef = self.obj.val

                antAxis = self.getAxesNames().index('ant')
                refSelection = self.selection[:]

                if refAnt == 'closest':
                    # put antenna axis first
                    dataVals = np.swapaxes(dataVals, 0, antAxis)

                    for i, antToRef in enumerate(self.getAxisValues('ant')):
                        # get the closest antenna
                        antDists = self.getSolset().getAntDist(antToRef) # this is a dict
                        for badAnt in self._getFullyFlaggedAnts(): del antDists[badAnt] # remove bad ants

                        refAnt = list(antDists.keys())[list(antDists.values()).index( sorted(antDists.values())[1] ) ] # get the second closest antenna (the first is itself)

                        refSelection[antAxis] = [self.getAxisValues('ant', ignoreSelection=True).tolist().index(refAnt)]
                        dataValsRef_i = self._applyAdvSelection(dataValsRef, refSelection)
                        dataValsRef_i = np.swapaxes(dataValsRef_i, 0, antAxis)
                        if weight:
                            dataVals[i][ dataValsRef_i[0] == 0. ] = 0.
                        else:
                            dataVals[i] -= dataValsRef_i[0]

                    dataVals = np.swapaxes(dataVals, 0, antAxis)
 
                else:
                    refSelection[antAxis] = [self.getAxisValues('ant', ignoreSelection=True).tolist().index(refAnt)]
                    dataValsRef = self._applyAdvSelection(dataValsRef, refSelection)
    
                    if weight:
                        dataVals[ np.repeat(dataValsRef, axis=antAxis, repeats=len(self.getAxisValues('ant'))) == 0. ] = 0.
                    else:
                        dataVals = dataVals - np.repeat(dataValsRef, axis=antAxis, repeats=len(self.getAxisValues('ant')))
                # if not weight and not self.getType() != 'tec' and not self.getType() != 'clock' and not self.getType() != 'tec3rd' and not self.getType() != 'rotationmeasure':
                #     dataVals = normalize_phase(dataVals)
        # CASE 2: Reference only to dir
        # TODO: should there be a warning if only direction is referenced but multipled ants are present?
        elif refDir and not refAnt:
            if not self.getType() in ['phase', 'scalarphase', 'rotation', 'tec', 'clock', 'tec3rd', 'rotationmeasure']:
                logging.error(
                    'Reference possible only for phase, scalarphase, clock, tec, tec3rd, rotation and rotationmeasure solution tables. Ignore referencing.')
            elif not 'dir' in self.getAxesNames():
                logging.error('Cannot find direction axis for referencing phases. Ignore referencing.')
            elif not refDir in self.getAxisValues('dir', ignoreSelection=True) and refDir != 'center':
                logging.error('Cannot find direction ' + refDir + '. Ignore referencing.')
            else:
                if self.useCache:
                    if weight:
                        dataValsRef = self.cacheWeight
                    else:
                        dataValsRef = self.cacheVal
                else:
                    if weight:
                        dataValsRef = self.obj.weight
                    else:
                        dataValsRef = self.obj.val

                dirAxis = self.getAxesNames().index('dir')
                refSelection = self.selection[:]

                if refDir == 'center':
                    # get the center (=closest to average) direction
                    dirsDict = self.getSolset().getSou()
                    meanDir = np.mean([dirsDict[k] for k in dirsDict], axis=0)
                    refDir, _ = min(dirsDict.items(), key=lambda kd: np.linalg.norm(kd[1] - meanDir))

                refSelection[dirAxis] = [self.getAxisValues('dir', ignoreSelection=True).tolist().index(refDir)]
                dataValsRef = self._applyAdvSelection(dataValsRef, refSelection)

                if weight:
                    dataVals[
                        np.repeat(dataValsRef, axis=dirAxis, repeats=len(self.getAxisValues('dir'))) == 0.] = 0.
                else:
                    dataVals = dataVals - np.repeat(dataValsRef, axis=dirAxis,
                                                    repeats=len(self.getAxisValues('dir')))
                # if not weight and not self.getType() != 'tec' and not self.getType() != 'clock' and not self.getType() != 'tec3rd' and not self.getType() != 'rotationmeasure':
                #     dataVals = normalize_phase(dataVals)

        # CASE 3: Reference to ant and to dir
        if refAnt and refDir:
            if not self.getType() in ['phase', 'scalarphase', 'rotation', 'tec', 'clock', 'tec3rd', 'rotationmeasure']:
                logging.error('Reference possible only for phase, scalarphase, clock, tec, tec3rd, rotation and rotationmeasure solution tables. Ignore referencing.')
            elif not 'ant' in self.getAxesNames():
                logging.error('Cannot find antenna axis for referencing phases. Ignore referencing.')
            elif not 'dir' in self.getAxesNames():
                logging.error('Cannot find direction axis for referencing phases. Ignore referencing.')
            elif not refAnt in self.getAxisValues('ant', ignoreSelection = True) and refAnt != 'closest':
                logging.error('Cannot find antenna '+refAnt+'. Ignore referencing.')
            elif refAnt == 'closest': # TODO: This needs to be implemented...
                logging.error('refAnt=\'closest\' is not supported (yet) when also referencing a direction.')
            elif not refDir in self.getAxisValues('dir', ignoreSelection=True) and refDir != 'center':
                logging.error('Cannot find direction ' + refDir + '. Ignore referencing.')
            else:
                if self.useCache:
                    if weight:
                        dataValsRef = self.cacheWeight
                    else:
                        dataValsRef = self.cacheVal
                else:
                    if weight:
                        dataValsRef = self.obj.weight
                    else:
                        dataValsRef = self.obj.val

                antAxis = self.getAxesNames().index('ant')
                dirAxis = self.getAxesNames().index('dir')
                refSelection = self.selection[:]

                if refDir == 'center':
                    # get the center (=closest to average) direction
                    dirsDict = self.getSolset().getSou()
                    meanDir = np.mean([dirsDict[k] for k in dirsDict], axis=0)
                    refDir, _ = min(dirsDict.items(), key=lambda kd: np.linalg.norm(kd[1] - meanDir))

                refSelection[antAxis] = [self.getAxisValues('ant', ignoreSelection=True).tolist().index(refAnt)]
                refSelection[dirAxis] = [self.getAxisValues('dir', ignoreSelection=True).tolist().index(refDir)]

                dataValsRef = self._applyAdvSelection(dataValsRef, refSelection)

                if weight:
                    dataValsRef = np.repeat(dataValsRef, dataVals.shape[antAxis], axis=antAxis)
                    dataValsRef = np.repeat(dataValsRef, dataVals.shape[dirAxis], axis=dirAxis)
                    dataVals[dataValsRef == 0.] = 0.

                else:
                    dataVals = dataVals - dataValsRef # np.expand_dims(dataValsRef, axis=(antAxis,dirAxis))

        if not retAxesVals:
            return dataVals

        axisVals = {}
        for axis in self.getAxesNames():
            axisVals[axis] = self.getAxisValues(axis)

        return dataVals, axisVals


    @deprecated_alias(reference='refAnt') # Add alias for backwards compatibility
    def getValuesIter(self, returnAxes=[], weight=False, refAnt=None, refDir=None):
        """
        Return an iterator which yields the values matrix (with axes = returnAxes) iterating along the other axes.
        E.g. if returnAxes are ['freq','time'], one gets a interetion over all the possible NxM
        matrix where N are the freq and M the time dimensions. The other axes are iterated in the getAxesNames() order.
        Note that all the data are fetched in memory before returning them one at a time. This is quicker.

        Parameters
        ----------
        returnAxes : list
            Axes of the returned array, all _others_ will be cycled on each element combinations.
        weight : bool, optional
            If true return also the weights, by default False.
        refAnt : str
            In case of phase solutions, reference to this station name.
        refDir : str
            In case of phase solutions, reference to this direction name.

        Returns
        -------
        1) data ndarray of dim=dim(returnAxes) and with the axes ordered as in getAxesNames()
        2) (if weight == True) weigth ndarray of dim=dim(returnAxes) and with the axes ordered as in getAxesNames()
        3) a dict with axis values in the form:
        {'axisname1':[axisvals1],'axisname2':[axisvals2],...}
        4) a selection which should be used to write this data back using a setValues()
        """
        if weight: weigthVals = self.getValues(retAxesVals=False, weight=True, refAnt=refAnt, refDir=refDir)
        dataVals = self.getValues(retAxesVals=False, weight=False, refAnt=refAnt, refDir=refDir)

        # get dimensions of non-returned axis (in correct order)
        iterAxesDim = [self.getAxisLen(axis) for axis in self.getAxesNames() if not axis in returnAxes]

        # generator to cycle over all the combinations of iterAxes
        # it "simply" gets the indexes of this particular combination of iterAxes
        # and use them to refine the selection.
        def g():
            for axisIdx in np.ndindex(tuple(iterAxesDim)):
                refSelection = []
                returnSelection = []
                thisAxesVals = {}
                i = 0
                for j, axisName in enumerate(self.getAxesNames()):
                    if axisName in returnAxes:
                        thisAxesVals[axisName] = self.getAxisValues(axisName)
                        # add a slice with all possible values (main selection is preapplied)
                        refSelection.append(slice(None))
                        # for the return selection use the "main" selection for the return axes
                        returnSelection.append(self.selection[j])
                    else:
                        #TODO: the iteration axes are not into a 1 element array, is it a problem?
                        thisAxesVals[axisName] = self.getAxisValues(axisName)[axisIdx[i]]
                        # add this index to the refined selection, this will return a single value for this axis
                        # an int is appended, this will remove an axis from the final data
                        refSelection.append(axisIdx[i])
                        # for the return selection use the complete axis and find the correct index
                        returnSelection.append( [self.getAxisValues(axisName, ignoreSelection=True).tolist().index(thisAxesVals[axisName])] )
                        i += 1

                # costly command
                data = dataVals[tuple(refSelection)]
                if weight:
                    weights = weigthVals[tuple(refSelection)]
                    yield (data, weights, thisAxesVals, returnSelection)
                else:
                    yield (data, thisAxesVals, returnSelection)

        return g()


    def addHistory(self, entry, date = True):
        """
        Adds entry to the table history with current date and time

        Since attributes cannot by default be larger than 64 kB, each
        history entry is stored in a separate attribute.

        Parameters
        ----------
        entry : str
            entry to add to history list
        """
        import datetime
        current_time = str(datetime.datetime.now()).split('.')[0]
        attrs = self.obj.val.attrs._f_list("user")
        nums = []
        for attr in attrs:
            try:
                if attr[:-3] == 'HISTORY':
                    nums.append(int(attr[-3:]))
            except:
                pass
        historyAttr = "HISTORY%03d" % min(list(set(range(1000)) - set(nums)))

        if date:
            entry = current_time + ": " + str(entry)
        else:
            entry = str(entry)

        self.obj.val.attrs[historyAttr] = entry.encode()


    def getHistory(self):
        """
        Get the soltab history.

        Returns
        -------
        str
            The table history as a string with each entry separated by newlines.
        """
        attrs = self.obj.val.attrs._f_list("user")
        attrs.sort()
        history_list = []
        for attr in attrs:
            if attr[:-3] == 'HISTORY':
                history_list.append(self.obj.val.attrs[attr].decode())
        if len(history_list) == 0:
            history_str = ""
        else:
            history_str = "\n".join(history_list)

        return history_str
