#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Retrieving and writing data in H5parm format

import os, sys, re
import numpy as np
import tables
import logging
import _version

# check for tables version
if int(tables.__version__.split('.')[0]) < 3:
    logging.critical('pyTables version must be >= 3.0.0, found: '+tables.__version__)
    sys.exit(1)

class h5parm():

    def __init__(self, h5parmFile, readonly = True, complevel = 5, complib='lzo'):
        """
        Keyword arguments:
        h5parmFile -- H5parm filename
        readonly -- if True the table is open in readonly mode (default=True)
        complevel -- compression level from 0 to 9 (default=5) when creating the file
        complib -- library for compression: lzo, zlib, bzip2 (default=lzo)
        """
        if os.path.isfile(h5parmFile):
            if tables.is_pytables_file(h5parmFile) == None:
                logging.critical('Wrong HDF5 format for '+h5parmFile+'.')
                raise Exception('Wrong HDF5 format for '+h5parmFile+'.')
            if readonly:
                logging.debug('Reading from '+h5parmFile+'.')
                self.H = tables.openFile(h5parmFile, 'r')
            else:
                logging.debug('Appending to '+h5parmFile+'.')
                self.H = tables.openFile(h5parmFile, 'r+')
        else:
            if readonly:
                raise Exception('Missing file '+h5parmFile+'.')
            else:
                logging.debug('Creating '+h5parmFile+'.')
                # add a compression filter
                f = tables.Filters(complevel=complevel, complib=complib)
                self.H = tables.openFile(h5parmFile, filters=f, mode='w')

        self.fileName = h5parmFile


    def __del__(self):
        """
        Flush and close the open table
        """
        self.H.close()


    def __str__(self):
        """
        Returns string with info about H5parm contents
        """
        return self.printInfo()


    def makeSolset(self, solsetName = None, addTables=True):
        """
        Create a new solset, if the provided name is not given or exists
        then it falls back on the first available sol###
        Keyword arguments:
        solset -- name of the solution set
        addTables -- if True (default) add antenna/direction/array tables
        """

        if type(solsetName) is str and not re.match(r'^[A-Za-z0-9_-]+$', solsetName):
            logging.warning('Solution-set '+solsetName+' contains unsuported characters. Use [A-Za-z0-9_-]. Switching to default.')
            solsetName = None

        if solsetName in self.getSolsets().keys():
            logging.warning('Solution-set '+solsetName+' already present. Switching to default.')
            solsetName = None

        if solsetName == None:
            solsetName = self._fisrtAvailSolsetName()

        logging.info('Creating a new solution-set: '+solsetName+'.')
        solset = self.H.create_group("/", solsetName)

        if addTables:
            # add antenna table
            logging.info('--Creating new antenna table.')
            descriptor = np.dtype([('name', np.str_, 16),('position', np.float32, 3)])
            soltab = self.H.createTable(solset, 'antenna', descriptor, \
                    title = 'Antenna names and positions', expectedrows = 40)
            soltab.attrs['h5parm_version'] = _version.__h5parmVersion__

            # add direction table
            logging.info('--Creating new source table.')
            descriptor = np.dtype([('name', np.str_, 16),('dir', np.float32, 2)])
            soltab = self.H.createTable(solset, 'source', descriptor, \
                    title = 'Source names and directions', expectedrows = 10)
            soltab.attrs['h5parm_version'] = _version.__h5parmVersion__

        return solset


    def getSolsets(self):
        """
        Return a dict with all the available solultion-sets (as Groups objects)
        """
        return self.H.root._v_groups


    def getSolset(self, solset = None):
        """
        Return a solultion-set as a Group object
        Keyword arguments:
        solset -- name of the solution set
        """
        if solset == None:
            raise Exception("Solution set not specified.")

        if not solset in self.getSolsets():
            logging.critical("Cannot find solset: "+solset+".")
            raise Exception("Cannot find solset: "+solset+".")

        return self.H.get_node('/',solset)


    def _fisrtAvailSolsetName(self):
        """
        Create and return the first available solset name which
        has the form of "sol###"
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
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        soltype -- solution-type (e.g. amplitude, phase)
        soltab -- the solution-table name (String) if not specified is generated from the solution-type
        axesNames -- list with the axes names
        axesVals -- list with the axes values
        chunkShape -- list with the chunk shape
        vals --
        weights --
        parmdbType -- original parmdb solution type
        """

        if soltype == None:
            raise Exception("Solution-type not specified while adding a solution-table.")

        # checks on the solset
        if solset == None:
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

        if soltabName == None:
            soltabName = self._fisrtAvailSoltabName(solset, soltype)

        logging.info('Creating a new solution-table: '+soltabName+'.')
        soltab = self.H.create_group("/"+solsetName, soltabName, title=soltype)
        soltab._v_attrs['parmdb_type'] = parmdbType

        # create axes
        assert len(axesNames) == len(axesVals)
        dim = []
        for i, axisName in enumerate(axesNames):
            axis = self.H.create_carray('/'+solsetName+'/'+soltabName, axisName,\
                    obj=axesVals[i], chunkshape=[len(axesVals[i])])
            axis.attrs['h5parm_version'] = _version.__h5parmVersion__
            dim.append(len(axesVals[i]))

        # check if the axes were in the proper order
        assert dim == list(vals.shape)
        assert dim == list(weights.shape)

        if chunkShape == None:
            chunkShape = 1+np.array(dim)/4

        # create the val/weight Carrays
        val = self.H.create_carray('/'+solsetName+'/'+soltabName, 'val', obj=vals, chunkshape=chunkShape)
        weight = self.H.create_carray('/'+solsetName+'/'+soltabName, 'weight', obj=weights, chunkshape=chunkShape)
        val.attrs['VERSION_H5PARM'] = _version.__h5parmVersion__
        val.attrs['AXES'] = ','.join([axisName for axisName in axesNames])
        weight.attrs['VERSION_H5PARM'] = _version.__h5parmVersion__
        weight.attrs['AXES'] = ','.join([axisName for axisName in axesNames])

        return soltab


    def getSoltabs(self, solset=None):
        """
        Return a dict {name1: object1, name2: object2, ...}
        of all the available solultion-tables into a specified solution-set
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        """
        if solset == None:
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
        if solset == None:
            raise Exception("Solution-set not specified while querying for solution-table.")
        if soltab == None:
            raise Exception("Solution-table not specified while querying for solution-table.")

        if not soltab in self.getSoltabs(solset):
            logging.critical("Solution-table "+soltab+" not found in solset "+solset+".")
            raise Exception("Solution-table "+soltab+" not found in solset "+solset+".")

        if type(solset) is str:
            solset = self.getSolset(solset)

        return self.H.get_node(solset, soltab)


    def _fisrtAvailSoltabName(self, solset=None, soltype=None):
        """
        Return the first available solset name which
        has the form of "soltypeName###"
        Keyword arguments:
        solset -- a solution-set name as Group instance
        soltype -- type of solution (amplitude, phase, RM, clock...) as a string
        """
        if solset == None:
            raise Exception("Solution-set not specified while querying for solution-tables list.")
        if soltype == None:
            raise Exception("Solution type not specified while querying for solution-tables list.")

        nums = []
        for soltab in self.getSoltabs(solset).keys():
            if re.match(r'^'+soltype+'[0-9][0-9][0-9]$', soltab):
                nums.append(int(soltab[-3:]))

        return soltype+"%03d" % min(list(set(range(1000)) - set(nums)))


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
        if solset == None:
            raise Exception("Solution-set not specified.")
        if type(solset) is str:
            solset = self.H.root._f_get_child(solset)

        sources = {}
        for x in solset.source:
            sources[x['name']] = x['dir']

        return sources


    def printInfo(self, filter=None):
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

            soltabs = self.getSoltabs(solset=solset_name)
            if len(soltabs) == 0:
                info += "\nNo tables\n"
            else:
                # For each table, print length of each axis and history of
                # operations applied to the table.
                for soltab_name in soltabs.keys():
                    try:
                        sf = solFetcher(soltabs[soltab_name])
                        axisNames = sf.getAxesNames()
                        axis_str_list = []
                        for axisName in axisNames:
                            nslots = len(sf.getAxisValues(axisName))
                            if nslots > 1:
                                pls = "s"
                            else:
                                pls = ""
                            axis_str_list.append("%i %s%s" % (nslots, axisName, pls))
                        info += "\nSolution table '%s': %s\n" % (soltab_name, ", ".join(axis_str_list))
                        history = sf.getHistory()
                        if history != "":
                            info += "\n" + 4*" " + "History:\n" + 4*" "
                            joinstr =  "\n" + 4*" "
                            info += joinstr.join(wrap(history)) + "\n"
                    except tables.exceptions.NoSuchNodeError:
                        info += "\nSolution table '%s': No valid data found\n" % (soltab_name)

        return info


class solHandler():
    """
    Generic class to principally handle selections
    Selections are:
    axisName = None # to select all
    axisName = xxx # to select ONLY that value for an axis
    axisName = [xxx, yyy, zzz] # to selct ONLY those values for an axis
    axisName = 'xxx' # regular expression selection
    axisName = {min: xxx} # to selct values grater or equal than xxx
    axisName = {max: yyy} # to selct values lower or equal than yyy
    axisName = {min: xxx, max: yyy} # to selct values greater or equal than xxx and lower or equal than yyy
    """
    def __init__(self, table, **args):
        """
        Keyword arguments:
        table -- table object
        **args -- used to create a selection
        """

        if not isinstance( table, tables.Group ):
            logging.error("Object must be initialized with a pyTables Table object.")
            return None

        self.t = table
        self.selection = {}
        self.setSelection(**args)


    def setSelection(self, **args):
        """
        set a default selection criteria.
        Keyword arguments:
        *args -- valid axes names of the form: pol='XX', ant=['CS001HBA','CS002HBA'], time={'min':1234,'max':'2345'}.
        """
        # create an initial selection which selects all values
        # any selection will modify only the slice relative to that axis
        self.selection = [slice(0,self.getAxisLen(axisName, ignoreSelection=True)) \
                                                for axisName in self.getAxesNames()]

        for axis, selVal in args.iteritems():
            if not axis in self.getAxesNames():
                logging.error("Cannot select on axis "+axis+", it doesn't exist.")
                return

            # if None continue and keep all the values
            if selVal == None: continue

            # find the index of the working axis
            idx = self.getAxesNames().index(axis)

            # string -> regular expression
            if type(selVal) is str:
                if self.getAxis(axis).atom.dtype.kind != 'S':
                    logging.error("Cannot select on axis "+axis+" with a regular expression.")
                    return
                self.selection[idx] = [i for i, item in enumerate(self.getAxisValues(axis)) if re.search(selVal, item)]

                # transform list of 1 element in a relative slice(), necessary when slicying and to always get an array back    
                if len(self.selection[idx]) == 1: self.selection[idx] = slice(self.selection[idx][0],self.selection[idx][0]+1,None)

            # dict -> min max
            elif type(selVal) is dict:
                axisVals = self.getAxisValues(axis)
                if 'min' in selVal and 'max' in selVal:
                    self.selection[idx] = slice(np.where(axisVals >= selVal['min'])[0][0],np.where(axisVals <= selVal['max'])[0][-1]+1)
                elif 'min' in selVal:
                    self.selection[idx] = slice(np.where(axisVals >= selVal['min'])[0][0],self.getAxisLen(axisName, ignoreSelection=True))
                elif 'max' in selVal:
                    self.selection[idx] = slice(0,np.where(axisVals <= selVal['max'])[0][-1]+1)
                else:
                    logging.error("Selection with a dict must have 'min' and/or 'max' entry.")
                    return

            # single val/list -> exact matching
            else:
                if not type(selVal) is list: selVal = [selVal]
                self.selection[idx] = [i for i, item in enumerate(self.getAxisValues(axis)) if item in selVal]
            
                # transform list of 1 element in a relative slice(), necessary when slicying and to always get an array back    
                if len(self.selection[idx]) == 1: self.selection[idx] = slice(self.selection[idx][0],self.selection[idx][0]+1,None)

            # if a selection return an empty list (maybe because of a wrong name), then use all values
            if type(self.selection[idx]) is list and len(self.selection[idx]) == 0:
                logging.error("Empty/wrong selection on axis \""+axis+"\". Use all available values.")
                self.selection[idx] = slice(0,self.getAxisLen(axisName, ignoreSelection=True))


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
        return self.t.val.attrs['AXES'].split(',')


    def getAxis(self, axis = None):
        """
        Return the axis istance for the corresponding name
        Keyword arguments:
        axis -- the name of the axis to be returned
        """
        if not axis in self.getAxesNames():
            logging.error("Cannot find "+axis+", it doesn't exist.")
            return None
        return self.t._f_get_child(axis)


    def getAxisLen(self, axis = None, ignoreSelection = False):
        """
        Return the axis lenght
        Keyword arguments:
        axis -- the name of the axis to be returned
        ignoreSelection -- if True returns the axis lenght without any selection active
        """
        if not axis in self.getAxesNames():
            logging.error("Cannot find "+axis+", it doesn't exist.")
            return None
        if ignoreSelection:
            return self.t._f_get_child(axis).nrows
        else:
            return len(self.getAxisValues(axis))


    def getAxisValues(self, axis='', ignoreSelection = False):
        """
        Return a list of all the possible values present along a specific axis (no duplicates)
        Keyword arguments:
        axis -- the axis name
        ignoreSelection -- if True returns the axis values without any selection active
        """

        if axis not in self.getAxesNames():
            logging.error('Axis \"'+axis+'\" not found.')
            return None

        axisIdx = self.getAxesNames().index(axis)
        if ignoreSelection:
            return np.copy(self.getAxis(axis)[:])
        else:
            return np.copy(self.getAxis(axis)[self.selection[axisIdx]])


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

    def __init__(self, table, **args):
        solHandler.__init__(self, table = table, **args)


    def setAxisValues(self, axis = None, vals = None):
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
        if weight: dataVals = self.t.weight
        else: dataVals = self.t.val

        dataVals[tuple(self.selection)] = vals


class solFetcher(solHandler):

    def __init__(self, table, **args):
        solHandler.__init__(self, table = table, **args)


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
            logging.error("Cannot find axis \""+axis+"\".")
            return None


    def getValues(self, retAxesVals = True, weight = False):
        """
        Creates a simple matrix of values. Fetching all selected rows into memory.
        Keyword arguments:
        retAxisVals -- if true returns also the axes vals as a dict of:
        {'axisname1':[axisvals1],'axisname2':[axisvals2],...}
        weight -- if true store in the weights instead that in the vals (default: False)
        """

        if weight: dataVals = self.t.weight
        else: dataVals = self.t.val

        # Use the slices set by setSelection to slice the data
        if not retAxesVals:
            return dataVals[tuple(self.selection)]

        axisVals = {}
        for axis in self.getAxesNames():
            axisVals[axis] = self.getAxisValues(axis)

        return dataVals[tuple(self.selection)], axisVals


    def getValuesIter(self, returnAxes= [], weight = False):
        """
        Return an iterator which yields the values matrix (with axes = returnAxes) iterating along the other axes.
        E.g. if returnAxes are ['freq','time'], one gets a interetion over all the possible NxM
        matrix where N are the freq and M the time dimensions. The other axes are iterated in the getAxesNames() order.
        Keyword arguments:
        returnAxes -- axes of the returned array, all others will be cycled
        weight -- if true store in the weights instead that in the vals (default: False)
        Return:
        1) ndarray of dim=dim(returnAxes) and with the axes ordered as in getAxesNames()
        2) a dict with axis values in the form:
        {'axisname1':[axisvals1],'axisname2':[axisvals2],...}
        """
        if weight: dataVals = self.t.weight[tuple(self.selection)]
        else: dataVals = self.t.val[tuple(self.selection)]

        # get dimensions of non-returned axis (in correct order)
        iterAxesDim = [self.getAxisLen(axis) for axis in self.getAxesNames() if not axis in returnAxes]

        # generator to cycle over all the combinations of iterAxes
        # it "simply" get the indexes of this particular combination of iterAxes
        # and use them to refine the selection.
        def g():
            for axisIdx in np.ndindex(tuple(iterAxesDim)):
                refSelection = []
                thisAxesVals = {}
                i = 0
                for axisName in self.getAxesNames():
                    if axisName in returnAxes:
                        thisAxesVals[axisName] = self.getAxisValues(axisName)
                        # add a slice with all possible values
                        refSelection.append(slice(0,self.getAxisLen(axisName),None))
                    else:
                        thisAxesVals[axisName] = self.getAxisValues(axisName)[axisIdx[i]]
                        # add this index to the refined selection, this will return a single value for this axis
                        refSelection.append(axisIdx[i])
                        i += 1
                # costly command
                data = dataVals[tuple(refSelection)]
                yield (data, thisAxesVals)

        return g()
