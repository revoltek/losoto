#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Retrieving and writing data in H5parm format

import os, sys, re, itertools
import numpy as np
import tables
import logging
import _version

# check for tables version
if int(tables.__version__.split('.')[0]) < 3:
    logging.critical('pyTables version must be >= 3.0.0, found: '+tables.__version__)
    sys.exit(1)

class h5parm( object ):

    def __init__(self, h5parmFile, readonly=True, complevel=0, complib='zlib'):
        """
        Create an h5parm object.

        Parameters
        ----------
        h5parmFile : str
            H5parm filename
        readonly : bool, optional
            if True the table is open in readonly mode, by default True
        complevel : int, optional
            compression level from 0 to 9 when creating the file, by default 5
        complib : str, optional
            library for compression: lzo, zlib, bzip2, by default zlib
        """

        self.H = None # variable to store the pytable object

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

        self.fileName = h5parmFile


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
        then it falls back on the first available sol###

        Parameters
        ----------
        solset : str
            name of the solution set
        addTables : bool, optional
            if True add antenna/direction/array tables, by default True

        Returns
        -------
        pytables Group
            newly created solset object
        """

        if type(solsetName) is str and not re.match(r'^[A-Za-z0-9_-]+$', solsetName):
            logging.warning('Solution-set '+solsetName+' contains unsuported characters. Use [A-Za-z0-9_-]. Switching to default.')
            solsetName = None

        if solsetName in self.getSolsets().keys():
            logging.warning('Solution-set '+solsetName+' already present. Switching to default.')
            solsetName = None

        if solsetName is None:
            solsetName = self._firstAvailSolsetName()

        logging.info('Creating a new solution-set: '+solsetName+'.')
        solset = self.H.create_group("/", solsetName)
        solset._f_setattr('h5parm_version', _version.__h5parmVersion__)

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

        return solset


    def getSolsets(self):
        """
        Returns
        -------
        dict
            a dict with all the available solultion-sets and relative Group object
        """
        return self.H.root._v_groups


    def getSolset(self, solset = None):
        """
        Parameters
        ----------
        solset : str
            name of the solution set

        Returns
        -------
        pytables Group
            return the solultion-set
        """
        if solset is None:
            raise Exception("Solution set not specified.")

        if not solset in self.getSolsets():
            logging.critical("Cannot find solset: "+solset+".")
            raise Exception("Cannot find solset: "+solset+".")

        return self.H.get_node('/',solset)


    def _firstAvailSolsetName(self):
        """
        Find the first available solset name which has the form of "sol###"

        Returns
        -------
        str
            solset name
        """
        nums = []
        for solset in self.getSolsets().keys():
            if re.match(r'^sol[0-9][0-9][0-9]$', solset):
                nums.append(int(solset[-3:]))

        return "sol%03d" % min(list(set(range(1000)) - set(nums)))


    def makeSoltab(self, solset=None, soltype=None, soltab=None,
            axesNames = [], axesVals = [], chunkShape=None, vals=None,
            weights=None, parmdbType=None):
        """
        Create a solution-table into a specified solution-set
        
        Parameters
        ----------
        solset : str or pytables Group
            solution-set name or a Group instance
        soltype : str
            solution-type (e.g. amplitude, phase)
        soltab : str, optional
            the solution-table name, if not specified is generated from the solution-type
        axesNames : list
            list with the axes names
        axesVals : list
            list with the axes values (each is a separate list)
        chunkShape : list, optional
            list with the chunk shape
        vals : numpy array
            array with shape given by the axesVals lenghts
        weights : numpy array
            same shape of the vals array
            0->FLAGGED, 1->MAX_WEIGHT
        parmdbType : str
            original parmdb solution type

        Returns
        -------
        pytables Group
            newly created soltab object
        """

        if soltype is None:
            raise Exception("Solution-type not specified while adding a solution-table.")

        # checks on the solset
        if solset is None:
            raise Exception("Solution-set not specified while adding a solution-table.")
        if type(solset) is str:
            solset = self.getSolset(solset)
        solsetName = solset._v_name

        if not solsetName in self.getSolsets().keys():
            raise Exception("Solution-set "+solsetName+" doesn't exist.")

        # checks on the soltab
        soltabName = soltab
        if type(soltabName) is str and not re.match(r'^[A-Za-z0-9_-]+$', soltabName):
            logging.warning('Solution-table '+soltabName+' contains unsuported characters. Use [A-Za-z0-9_-]. Switching to default.')
            soltabName = None

        if soltabName in self.getSoltabs(solset).keys():
            logging.warning('Solution-table '+soltabName+' already present. Switching to default.')
            soltabName = None

        if soltabName is None:
            soltabName = self._fisrtAvailSoltabName(solset, soltype)

        logging.info('Creating a new solution-table: '+soltabName+'.')
        soltab = self.H.create_group("/"+solsetName, soltabName, title=soltype)
        soltab._v_attrs['parmdb_type'] = parmdbType

        # create axes
        assert len(axesNames) == len(axesVals)
        dim = []

#        newChunkShape = []
        for i, axisName in enumerate(axesNames):
            #axis = self.H.create_carray('/'+solsetName+'/'+soltabName, axisName,\
            #        obj=axesVals[i], chunkshape=[len(axesVals[i])])
            axis = self.H.create_array('/'+solsetName+'/'+soltabName, axisName, obj=axesVals[i])
            dim.append(len(axesVals[i]))

        # check if the axes were in the proper order
        assert dim == list(vals.shape)
        assert dim == list(weights.shape)

        # create the val/weight Carrays
        #val = self.H.create_carray('/'+solsetName+'/'+soltabName, 'val', obj=vals.astype(np.float64), chunkshape=None, atom=tables.Float64Atom())
        #weight = self.H.create_carray('/'+solsetName+'/'+soltabName, 'weight', obj=weights.astype(np.float16), chunkshape=None, atom=tables.Float16Atom())
        # array do not have compression but are much faster
        val = self.H.create_array('/'+solsetName+'/'+soltabName, 'val', obj=vals.astype(np.float64), atom=tables.Float64Atom())
        weight = self.H.create_array('/'+solsetName+'/'+soltabName, 'weight', obj=weights.astype(np.float16), atom=tables.Float16Atom())
        val.attrs['AXES'] = ','.join([axisName for axisName in axesNames])
        weight.attrs['AXES'] = ','.join([axisName for axisName in axesNames])

        return soltab
    

    def _fisrtAvailSoltabName(self, solset=None, soltype=None):
        """
        Return the first available solset name which
        has the form of "soltypeName###"
        Keyword arguments:
        solset -- a solution-set name as Group instance
        soltype -- type of solution (amplitude, phase, RM, clock...) as a string
        """
        if solset is None:
            raise Exception("Solution-set not specified while querying for solution-tables list.")
        if soltype is None:
            raise Exception("Solution type not specified while querying for solution-tables list.")

        nums = []
        for soltab in self.getSoltabs(solset).keys():
            if re.match(r'^'+soltype+'[0-9][0-9][0-9]$', soltab):
                nums.append(int(soltab[-3:]))

        return soltype+"%03d" % min(list(set(range(1000)) - set(nums)))


    def delSoltab(self, solset=None, soltab=None):
        """
        Delete a solution-table of a specific solution-set
        Keyword arguments:
        solset -- a solution-set name (String) or instance (required if soltab is a string)
        soltab -- a solution-table name (String) or instance
        """
        if soltab is None:
            raise Exception("Solution-table not specified while deleting a solution-table.")

        if type(soltab) is str:
            soltabobj = self.getSoltab(solset, soltab)

        soltabobj._f_remove(recursive=True, force=True)
        logging.info("Soltab \""+soltab+"\" deleted.")


    def getSoltabs(self, solset=None):
        """
        Return a dict {name1: object1, name2: object2, ...}
        of all the available solultion-tables into a specified solution-set
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        """
        if solset is None:
            raise Exception("Solution-set not specified while querying for solution-tables list.")
        if type(solset) is str:
            solset = self.getSolset(solset)

        return solset._v_groups


    def getSoltab(self, solset=None, soltab=None):
        """
        Return a specific solution-table (object) of a specific solution-set
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        soltab -- a solution-table name (String)
        """
        if solset is None:
            raise Exception("Solution-set not specified while querying for solution-table.")
        if soltab is None:
            raise Exception("Solution-table not specified while querying for solution-table.")

        if not soltab in self.getSoltabs(solset):
            logging.critical("Solution-table "+soltab+" not found in solset "+solset+".")
            raise Exception("Solution-table "+soltab+" not found in solset "+solset+".")

        if type(solset) is str:
            solset = self.getSolset(solset)

        return self.H.get_node(solset, soltab)


    def getAnt(self, solset):
        """
        Return a dict of all available antennas
        in the form {name1:[position coords],name2:[position coords],...}
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        """
        if solset == None:
            raise Exception("Solution-set not specified.")
        if type(solset) is str:
            solset = self.H.root._f_get_child(solset)

        ants = {}
        for x in solset.antenna:
            ants[x['name']] = x['position']

        return ants


    def getSou(self, solset):
        """
        Return a dict of all available sources
        in the form {name1:[ra,dec],name2:[ra,dec],...}
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        """
        if solset is None:
            raise Exception("Solution-set not specified.")
        if type(solset) is str:
            solset = self.H.root._f_get_child(solset)

        sources = {}
        for x in solset.source:
            sources[x['name']] = x['dir']

        return sources


    def printInfo(self, filter=None, verbose=False):
        """
        Returns string with info about H5parm contents
        """
        from itertools import izip_longest

        def grouper(n, iterable, fillvalue=' '):
            """
            Groups iterables into specified groups

            Keyword arguments:
            n -- number of iterables to group
            iterable -- iterable to group
            fillvalue -- value to use when to fill blanks in output groups

            Example:
            grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
            """
            args = [iter(iterable)] * n
            return izip_longest(fillvalue=fillvalue, *args)

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
            for solset_name in solsets.keys():
                if not re.search(filter, solset_name):
                    keys_to_remove.append(solset_name)
            for key in keys_to_remove:
                solsets.pop(key)

        if len(solsets) == 0:
            info += "\nNo solution sets found.\n"
            return info
        solset_names = solsets.keys()
        solset_names.sort()
        # delete axes value file if already present
        if verbose and os.path.exists(self.fileName+'-axes_values.txt'):
                logging.warning('Overwriting '+self.fileName+'-axes_values.txt')
                os.system('rm '+self.fileName+'-axes_values.txt')
        # For each solution set, list solution tables, sources, and antennas
        for solset_name in solset_names:
            info += "\nSolution set '%s':\n" % solset_name
            info += "=" * len(solset_name) + "=" * 16 + "\n\n"

            # Print direction (source) names
            sources = self.getSou(solset_name).keys()
            sources.sort()
            info += "Directions: "
            for src_name in sources:
                info += "%s\n            " % src_name

            # Print station names
            antennas = self.getAnt(solset_name).keys()
            antennas.sort()
            info += "\nStations: "
            for ant1, ant2, ant3, ant4 in grouper(4, antennas):
                info += "{0:<10s} {1:<10s} {2:<10s} {3:<10s}\n          ".format(ant1, ant2, ant3, ant4)

            # For each table, print length of each axis and history of
            # operations applied to the table.
            if verbose:
                logging.warning('Axes values saved in '+self.fileName+'-axes_values.txt')
                f = file(self.fileName+'-axes_values.txt','a')
            for soltab_name, soltab in self.getSoltabs(solset=solset_name).iteritems():
                try:
                    if verbose:
                        f.write("### /"+solset_name+"/"+soltab_name+"\n")
                    logging.debug('Fetching info for '+soltab_name+'.')
                    sf = solFetcher(soltab)
                    axisNames = sf.getAxesNames()
                    #print self.getSoltabs(solset=solset_name)
                    axis_str_list = []
                    for axisName in axisNames:
                        info
                        nslots = sf.getAxisLen(axisName)
                        if nslots > 1:
                            pls = "s"
                        else:
                            pls = ""
                        axis_str_list.append("%i %s%s" % (nslots, axisName, pls))
                        if verbose:
                            f.write(axisName+": ")
                            vals = sf.getAxisValues(axisName)
                            # ugly hardcoded workaround to print all the important decimal values for time/freq
                            if axisName == 'freq': f.write(" ".join(["{0:.8f}".format(v) for v in vals])+"\n\n")
                            elif axisName == 'time': f.write(" ".join(["{0:.7f}".format(v) for v in vals])+"\n\n")
                            else: f.write(" ".join(["{}".format(v) for v in vals])+"\n\n")
                    info += "\nSolution table '%s' (type: %s): %s\n" % (soltab_name, sf.getType(), ", ".join(axis_str_list))
                    weights = sf.getValues(weight = True, retAxesVals = False)
                    info += 'Flagged data %.3f%%\n' % (100.*np.sum(weights==0)/len(weights.flat))
                    history = sf.getHistory()
                    if history != "":
                        info += "\n" + 4*" " + "History:\n" + 4*" "
                        joinstr =  "\n" + 4*" "
                        info += joinstr.join(wrap(history)) + "\n"
                    del sf
                except tables.exceptions.NoSuchNodeError:
                    info += "\nSolution table '%s': No valid data found\n" % (soltab_name)

            if verbose:
                    f.close()
        return info


class solHandler( object ):
    """
    Generic #class to principally handle selections
    Selections are:
    axisName = None # to select all
    axisName = xxx # to select ONLY that value for an axis
    axisName = [xxx, yyy, zzz] # to selct ONLY those values for an axis
    axisName = 'xxx' # regular expression selection
    axisName = {min: xxx} # to selct values grater or equal than xxx
    axisName = {max: yyy} # to selct values lower or equal than yyy
    axisName = {min: xxx, max: yyy} # to selct values greater or equal than xxx and lower or equal than yyy
    """
    def __init__(self, table, useCache = False, **args):
        """
        Keyword arguments:
        table -- table object or solset/soltab string
        useCache -- cache all data in memory
        **args -- used to create a selection
        """

        if not isinstance( table, tables.Group ):
            logging.error("Object must be initialized with a pyTables Table object.")
            sys.exit(1)

        self.t = table
        # set axes names once to speed up calls
        self.axesNames = table.val.attrs['AXES'].split(',')
        # set axes values once to speed up calls (a bit of memory usage though)
        self.axes = {}
        for axis in self.getAxesNames():
            self.axes[axis] = table._f_get_child(axis)

        self.selection = {}
        self.setSelection(**args)

        self.useCache = useCache
        if self.useCache:
            logging.debug("Caching...")
            self.setCache(np.copy(self.t.val), np.copy(self.t.weight))


    def setCache(self, val, weight):
        """
        Set cache value
        """
        self.cacheVal = val
        self.cacheWeight = weight


    def getAddress(self):
        """
        return the solset/soltab address of self.t as a string
        """
        return self.t._v_pathname[1:]


    def setSelection(self, **args):
        """
        set a default selection criteria.
        Keyword arguments:
        *args -- valid axes names of the form: pol='XX', ant=['CS001HBA','CS002HBA'], time={'min':1234,'max':'2345','step':4}.
        """
        # create an initial selection which selects all values
        # any selection will modify only the slice relative to that axis
        self.selection = [slice(0,self.getAxisLen(axisName, ignoreSelection=True)) \
                                                for axisName in self.getAxesNames()]

        for axis, selVal in args.iteritems():
            # if None continue and keep all the values
            if selVal is None: continue

            if not axis in self.getAxesNames():
                logging.error("Cannot select on axis "+axis+", it doesn't exist. Ignored.")
                continue

            # find the index of the working axis
            idx = self.getAxesNames().index(axis)

            # string -> regular expression
            if type(selVal) is str:
                if not self.getAxisType(axis).char is 'S':
                    logging.warning("Cannot select on axis \""+axis+"\" with a regular expression. Use all available values.")
                    continue
                self.selection[idx] = [i for i, item in enumerate(self.getAxisValues(axis)) if re.search(selVal, item)]

                # transform list of 1 element in a relative slice(), necessary when slicying and to always get an array back
                if len(self.selection[idx]) == 1: self.selection[idx] = slice(self.selection[idx][0],self.selection[idx][0]+1)

            # dict -> min max
            elif type(selVal) is dict:
                axisVals = self.getAxisValues(axis)
                if 'min' in selVal and selVal['min'] > np.max(axisVals):
                    logging.error("Selection with min > than maximum value. Use all available values.")
                    continue
                if 'max' in selVal and selVal['max'] < np.min(axisVals):
                    logging.error("Selection with max < than minimum value. Use all available values.")
                    continue
                if 'min' in selVal and 'max' in selVal:
                    self.selection[idx] = slice(np.where(axisVals >= selVal['min'])[0][0], np.where(axisVals <= selVal['max'])[0][-1]+1)
                elif 'min' in selVal:
                    self.selection[idx] = slice(np.where(axisVals >= selVal['min'])[0][0], None)
                elif 'max' in selVal:
                    self.selection[idx] = slice(0, np.where(axisVals <= selVal['max'])[0][-1]+1)
                else:
                    logging.error("Selection with a dict must have 'min' and/or 'max' entry. Use all available values.")
                    continue
                if 'step' in selVal:
                    self.selection[idx] = slice(self.selection[idx].start, self.selection[idx].stop, selVal['step'])

            # single val/list -> exact matching
            else:
                if type(selVal) is np.array or type(selVal) is np.ndarray: selVal = selVal.tolist()
                if not type(selVal) is list: selVal = [selVal]
                # convert to correct data type (from parset everything is a str)
                selVal = np.array(selVal, dtype=self.getAxisType(axis))

                if len(selVal) == 1:
                    # speedup in the common case of a single value
                    self.selection[idx] = [self.getAxisValues(axis).tolist().index(selVal)]
                else:
                    self.selection[idx] = [i for i, item in enumerate(self.getAxisValues(axis)) if item in selVal]

                # transform list of 1 element in a relative slice(), necessary when slicying and to always get an array back
                if len(self.selection[idx]) == 1: self.selection[idx] = slice(self.selection[idx][0], self.selection[idx][0]+1)
                # transform list of continuous numbers in slices (faster)
                elif len(self.selection[idx]) != 0 and len(self.selection[idx])-1 == self.selection[idx][-1] - self.selection[idx][0]:
                    self.selection[idx] = slice(self.selection[idx][0], self.selection[idx][-1]+1)

            # if a selection return an empty list (maybe because of a wrong name), then use all values
            if type(self.selection[idx]) is list and len(self.selection[idx]) == 0:
                logging.warning("Empty/wrong selection on axis \""+axis+"\". Use all available values.")
                self.selection[idx] = slice(None)


    def getType(self):
        """
        return the type of the solution-tables (it is stored in an attrs)
        """

        return self.t._v_title


    def getAxesNames(self):
        """
        Return a list with all the axis names in the correct order for
        slicing the getValuesGrid() reurned matrix.
        """

        #return self.t.val.attrs['AXES'].split(',') # slower
        return self.axesNames[:]


    def getAxis(self, axis=None):
        """
        Return the axis istance for the corresponding name
        Keyword arguments:
        axis -- the name of the axis to be returned
        """
        try:
            return self.axes[axis]
        except:
            logging.error("Cannot find "+axis+", it doesn't exist.")
            return None


    def getAxisLen(self, axis=None, ignoreSelection=False):
        """
        Return the axis lenght
        Keyword arguments:
        axis -- the name of the axis to be returned
        ignoreSelection -- if True returns the axis lenght without any selection active
        """
        return len(self.getAxisValues(axis, ignoreSelection = ignoreSelection))


    def getAxisType(self, axis=None):
        """
        Return the axis dtype
        Keyword arguments:
        axis -- the name of the axis whose type is returned
        """
        try:
            return self.t._f_get_child(axis).dtype
        except:
            logging.error('Axis \"'+axis+'\" not found.')
            return None


    def getAxisValues(self, axis='', ignoreSelection=False):
        """
        Return a copy of all the possible values present along a specific axis (no duplicates)
        Keyword arguments:
        axis -- the axis name
        ignoreSelection -- if True returns the axis values without any selection active
        """
        try:
            if ignoreSelection:
                return self.getAxis(axis)[:]
            else:
                axisIdx = self.getAxesNames().index(axis)
                return self.getAxis(axis)[self.selection[axisIdx]]
        except:
            logging.error('Axis \"'+axis+'\" not found.')
            return None


    def addHistory(self, entry=""):
        """
        Adds entry to the table history with current date and time

        Since attributes cannot by default be larger than 64 kB, each
        history entry is stored in a separate attribute.

        Keyword arguments:
        entry -- string to add to history list
        """
        import datetime
        current_time = str(datetime.datetime.now()).split('.')[0]
        attrs = self.t.val.attrs._f_list("user")
        nums = []
        for attr in attrs:
            try:
                if attr[:-3] == 'HISTORY':
                    nums.append(int(attr[-3:]))
            except:
                pass
        historyAttr = "HISTORY%03d" % min(list(set(range(1000)) - set(nums)))

        self.t.val.attrs[historyAttr] = current_time + ": " + str(entry)


    def getHistory(self):
        """
        Returns the table history as a string with each entry separated by
        newlines
        """
        attrs = self.t.val.attrs._f_list("user")
        attrs.sort()
        history_list = []
        for attr in attrs:
            if attr[:-3] == 'HISTORY':
                history_list.append(self.t.val.attrs[attr])
        if len(history_list) == 0:
            history_str = ""
        else:
            history_str = "\n".join(history_list)

        return history_str


class solWriter(solHandler):

    def __init__(self, table, useCache = False, **args):
        """
        useCache -- write the data on a local copy of the table,
        use flush() to write them on the disk, speeds up writing
        """
        solHandler.__init__(self, table=table, useCache=useCache, **args)

    def setAxisValues(self, axis=None, vals=None):
        """
        Set the value of a specific axis
        Keyword arguments:
        axis -- the axis name
        vals -- the values
        """

        if axis not in self.getAxesNames():
            logging.error('Axis \"'+axis+'\" not found.')

        axisIdx = self.getAxesNames().index(axis)
        self.getAxis(axis)[self.selection[axisIdx]] = vals


    def setValues(self, vals, weight = False):
        """
        Save values in the val grid
        Keyword arguments:
        vals -- values to write as an n-dimentional array which match the selection dimention
        weight -- if true store in the weights instead that in the vals (default: False)
        """
        if self.useCache:
            if weight: dataVals = self.cacheWeight
            else: dataVals = self.cacheVal
        else:
            if weight: dataVals = self.t.weight
            else: dataVals = self.t.val

        # NOTE: pytables has a nasty limitation that only one list can be applied when selecting.
        # Conversely, one can apply how many slices he wants.
        # Single values/contigous values are converted in slices in h5parm.
        # This try/except implements a workaround for this limitation. Once the pytables will be updated, the except can be removed.
        try:
            # the float check allows quick reset of large arrays to a single value
            if isinstance(vals, (np.floating, float)): dataVals[tuple(self.selection)] = vals
            # the reshape is needed when saving e.g. [512] (vals shape) into [512,1,1] (selection output)
            else: dataVals[tuple(self.selection)] = np.reshape(vals, dataVals[tuple(self.selection)].shape)
        except:
            logging.debug('Optimizing selection writing '+str(self.selection))
            selectionListsIdx = [i for i, s in enumerate(self.selection) if type(s) is list]
            subSelection = self.selection[:]
            # create a subSelection also for the "vals" array
            subSelectionForVals = [slice(None) for i in xrange(len(subSelection))]
            # cycle across lists and save data index by index
            for selectionListValsIter in itertools.product(*[self.selection[selectionListIdx] for selectionListIdx in selectionListsIdx[1:]]):
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
        self.t.weight[:] = self.cacheWeight
        self.t.val[:] = self.cacheVal


class solFetcher(solHandler):

    def __init__(self, table, useCache = False, **args):
        solHandler.__init__(self, table = table, useCache = useCache, **args)

    def __getattr__(self, axis):
        """
        links any attribute with an "axis name" to getValuesAxis("axis name")
        also links val and weight to the relative arrays
        Keyword arguments:
        axis -- the axis name
        """
        if not axis in self.getAxesNames() and (axis != 'val' and axis != 'weight'):
            logging.error('Axis \"'+axis+'\" not found.')
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


    def getValues(self, retAxesVals=True, weight=False, reference=None, referencePol=None):
        """
        Creates a simple matrix of values. Fetching a copy of all selected rows into memory.
        Keyword arguments:
        retAxisVals -- if true returns also the axes vals as a dict of:
        {'axisname1':[axisvals1],'axisname2':[axisvals2],...}
        weight -- if true get the weights instead that the vals (default: False)
        reference -- in case of phase solutions, reference to this station name
        referencePol -- reference to selected pol of reference station
        Return:
        A numpy ndarrey (values or weights depending on parameters)
        If selected, returns also the axes values
        """

        if self.useCache:
            if weight: dataVals = self.cacheWeight
            else: dataVals = self.cacheVal
        else:
            if weight: dataVals = self.t.weight
            else: dataVals = self.t.val

        # apply the self.selection
        # NOTE: pytables has a nasty limitation that only one list can be applied when selecting.
        # Conversely, one can apply how many slices he wants.
        # Single values/contigous values are converted in slices in h5parm.
        # This try/except implements a workaround for this limitation. Once the pytables will be updated, the except can be removed.
        try:
            dataVals = dataVals[tuple(self.selection)]
        except:
            logging.debug('Optimizing selection reading '+str(self.selection))
            # for performances is important to minimize the fetched data
            # convert all lists other then the first (pytables allowes one!) to slices
            selectionListsIdx = [i for i, s in enumerate(self.selection) if type(s) is list]
            firstSelection = self.selection[:]
            for i in selectionListsIdx[1:]:
                firstSelection[i] = slice(None)
            # create a second selection using np.ix_
            secondSelection = []
            for i, sel in enumerate(self.selection):
                if i == selectionListsIdx[0]: secondSelection.append(xrange(self.getAxisLen(self.getAxesNames()[i], ignoreSelection=False)))
                elif type(sel) is list: secondSelection.append(sel)
                elif type(sel) is slice: secondSelection.append(xrange(self.getAxisLen(self.getAxesNames()[i], ignoreSelection=False)))
            dataVals = dataVals[tuple(firstSelection)][np.ix_(*secondSelection)]

        if reference is not None:
            # TODO: flag when reference is flagged?
            if self.getType() != 'phase' and self.getType() != 'scalarphase' and self.getType() != 'rotation' and self.getType() != 'tec' and self.getType() != 'clock' and self.getType() != 'tec3rd':
                logging.error('Reference possible only for phase, scalarphase, clock, tec, tec3rd, and rotation solution tables. Ignore referencing.')
            elif not 'ant' in self.getAxesNames():
                logging.error('Cannot find antenna axis for referencing phases. Ignore referencing.')
            elif not reference in self.getAxisValues('ant', ignoreSelection = True):
                logging.error('Cannot find antenna '+reference+'. Ignore referencing.')
            else:
                selection_stored = np.copy(self.selection)
                antAxis = self.getAxesNames().index('ant')
                self.selection[antAxis] = [list(self.getAxisValues('ant', ignoreSelection=True)).index(reference)]
                dataValsRef = self.getValues(retAxesVals=False, reference=None)
                if referencePol is not None:
                    polAxis = self.getAxesNames().index('pol')
                    # put pol axis at the beginning
                    dataValsRef = np.swapaxes(dataValsRef,0,polAxis)
                    # find reference pol index
                    polValIdx = list(self.getAxisValues('pol')).index(referencePol)
                    # set all polarisations equal to the ref pol, so at subtraction everything will be referenced only to that pol
                    for i in xrange(len(self.getAxisValues('pol'))):
                        dataValsRef[i] = dataValsRef[polValIdx]
                    # put pol axis back in place
                    dataValsRef = np.swapaxes(dataValsRef,0,polAxis)
                self.selection = selection_stored
                dataVals = dataVals - np.repeat(dataValsRef, axis=antAxis, repeats=len(self.getAxisValues('ant')))
                if not self.getType() != 'tec' and not self.getType() != 'clock':
                    # Convert to range [-2*pi, 2*pi].
                    dataVals = np.fmod(dataVals, 2.0 * np.pi)
                    # Convert to range [-pi, pi]
                    dataVals[dataVals < -np.pi] += 2.0 * np.pi
                    dataVals[dataVals > np.pi] -= 2.0 * np.pi

        if not retAxesVals:
            return dataVals

        axisVals = {}
        for axis in self.getAxesNames():
            axisVals[axis] = self.getAxisValues(axis)

        return dataVals, axisVals


    def getValuesIter(self, returnAxes=[], weight=False, reference=None, referencePol=None):
        """
        Return an iterator which yields the values matrix (with axes = returnAxes) iterating along the other axes.
        E.g. if returnAxes are ['freq','time'], one gets a interetion over all the possible NxM
        matrix where N are the freq and M the time dimensions. The other axes are iterated in the getAxesNames() order.
        Note that all the data are fetched in memory before returning them one at a time. This is quick.
        Keyword arguments:
        returnAxes -- axes of the returned array, all others will be cycled

        weight -- if true return also the weights (default: False)
        Return:
        1) data ndarray of dim=dim(returnAxes) and with the axes ordered as in getAxesNames()
        2) (if weight == True) weigth ndarray of dim=dim(returnAxes) and with the axes ordered as in getAxesNames()
        3) a dict with axis values in the form:
        {'axisname1':[axisvals1],'axisname2':[axisvals2],...}
        4) a selection which should be used to write this data back using a solWriter
        """
        if weight: weigthVals = self.getValues(retAxesVals=False, weight=True, reference=None)
        dataVals = self.getValues(retAxesVals=False, weight=False, reference=reference, referencePol=referencePol)

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
                        # add a slice with all possible values (main selection will be preapplied)
                        refSelection.append(slice(0,self.getAxisLen(axisName), None))
                        # for the return selection use the "main" selection for the return axes
                        returnSelection.append(self.selection[j])
                    else:
                        #TODO: the iteration axes are not into a 1 element array, is it a problem?
                        thisAxesVals[axisName] = self.getAxisValues(axisName)[axisIdx[i]]
                        # add this index to the refined selection, this will return a single value for this axis
                        refSelection.append(axisIdx[i])
                        # for the return selection use the complete axis and find the correct index
                        returnSelection.append( list(np.where( self.getAxisValues(axisName, ignoreSelection=True) == thisAxesVals[axisName] )[0]) )
                        i += 1
                # costly command
                data = dataVals[tuple(refSelection)]
                if weight:
                    weights = weigthVals[tuple(refSelection)]
                    yield (data, weights, thisAxesVals, returnSelection)
                else:
                    yield (data, thisAxesVals, returnSelection)

        return g()
