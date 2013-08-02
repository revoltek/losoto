#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Retrieving and writing data in H5parm format

import os, sys
import tables
import logging
import losoto._version

# check for tables version
if int(tables.__version__.split('.')[0]) < 3:
    logging.critical('pyTables version must be >= 3.0.0, found: '+tables.__version__)
    sys.exit(1)

class h5parm():

    def __init__(self, h5parmFile, readonly = True, complevel = 9):
        """
        Keyword arguments:
        h5parmFile -- H5parm filename
        readonly -- if True the table is open in readonly mode (default=True)
        complevel -- compression level from 0 to 9 (default=9) when creating the file
        """
        if os.path.isfile(h5parmFile):
            if readonly:
                self.H = tables.openFile(h5parmFile, 'r')
            else:
                self.H = tables.openFile(h5parmFile, 'a')
        else:
            if readonly:
                raise Exception('Missing file '+h5parmFile+'.')
            else:
                # add a compression filter
                f = tables.Filters(complevel=complevel, complib='zlib')
                self.H = tables.openFile(h5parmFile, filters=f, mode='w')
        
        # if the file is new add the version of the h5parm
        # in losoto._version.__h5parmVersion__


    def __del__(self):
        """
        Flush and close the open table
        """
        self.H.close()


    def makeSolset(self, solsetName = ''):
        """
        Create a new solset, if the provided name is not given or exists
        then it falls back on the first available sol###
        """

        if solsetName in self.getSolsets().keys():
            logging.error('Solution set '+solsetName+' already present.')
            solsetName = ''


        if solsetName == '':
            solsetName = self._fisrtAvailSolsetName()
        
        logging.info('Creating new solution-set '+solsetName+'.')
        return self.H.create_group("/", solsetName)


    def getSolsets(self):
        """
        Return a dict with all the available solultion-sets (as a _ChildrenDict)
        """
        return self.H.root._v_children


    def _fisrtAvailSolsetName(self):
        """
        Create and return the first available solset name which
        has the form of "sol###"
        """
        nums = []
        for solset in self.getSolsets().keys():
            try:
                if solset[0:3] == 'sol':
                    nums.append(int(solset[3:6]))
            except:
                pass

        return "sol%03d" % min(list(set(range(1000)) - set(nums)))


    def makeSoltab(self, solset=None, soltype=None, descriptor={}):
        """
        Create a solution-table into a specified solution-set
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        soltype -- solution type (e.g. amplitude, phase)
        """
        if solset == None:
            raise Exception("Solution set not specified while adding a solution-table.")
        if soltype == None:
            raise Exception("Solution type not specified while adding a solution-table.")
        
        if type(solset) is str:
            solset = self.H.root._f_get_child(solset)

        soltabName = self._fisrtAvailSoltabName(solset, soltype)
        logging.info('Creating new solution-table '+soltabName+'.')

        return self.H.createTable(solset, soltabName, descriptor, soltype)


    def getSoltabs(self, solset=None):
        """
        Return a dict {name1: object1, name2: object2, ...}
        of all the available solultion-tables into a specified solution-set
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        Output: 
        A dict of all available solultion-tables 
        """
        if solset == None:
            raise Exception("Solution set not specified while querying for solution-tables list.")
        if type(solset) is str:
            solset = self.H.root._f_get_child(solset)

        soltabs = {}
        for soltabName, soltab in solset._v_children.iteritems():
            if not (soltabName == 'antenna' or soltabName == 'source'):
                soltabs[soltabName] = soltab

        return soltabs


    def getSoltab(self, solset=None, soltab=None):
        """
        Return a specific solution-table of a specific solution-set
        Keyword arguments:
        solset -- a solution-set name (String) or a Group instance
        soltab -- a solution-table name (String)
        """
        if solset == None:
            raise Exception("Solution-set not specified.")
        if soltab == None:
            raise Exception("Solution-table not specified.")

        if type(solset) is str:
            solset = self.H.root._f_get_child(solset)

        return solset._f_get_child(soltab)


    def _fisrtAvailSoltabName(self, solset=None, soltype=None):
        """
        Create and return the first available solset name which
        has the form of "sol###"
        Keyword arguments:
        solset -- a solution-set name as Group instance
        soltype -- type of solution (amplitude, phase, RM, clock...) as a string
        """
        if solset == None:
            raise Exception("Solution-set not specified while querying for solution-tables list.")
        if soltype == None:
            raise Exception("Solution type not specified while querying for solution-tables list.")

        nums = []
        for soltab in self.getSoltabs(solset):
            try:
                if soltab[-4:] == soltype:
                    nums.append(int(soltab[-4:]))
            except:
                pass

        return soltype+"%03d" % min(list(set(range(1000)) - set(nums)))


    def addRow(self, soltab=None, val=[]):
        """
        Add a single row to the given soltab
        Keyword arguments:
        soltab -- a solution-table instance
        val -- a list of all the field to insert, the order is important!
        """
        if soltab == None:
            raise Exception("Solution-table not specified while adding a new row.")

        soltab.append(val)


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


class solFetcher():

    def __init__(self, table, selection = ''):
        """
        Keyword arguments:
        tab -- table object
        selection -- a selection on the axis of the type "(ant == 'CS001LBA') & (pol == 'XX')"
        """
        
        self.t = table
        self.selection = selection

    def setSelection(self, selection = ''):
        """
        set a default selection criteria.
        Keyword arguments:
        selection -- a selection on the axis of the type "(ant == 'CS001LBA') & (pol == 'XX')"
        """
        self.selection = selection


    def getType(self):
        """
        return the type of the solution-tables which is stored in the title
        """

        return self.t._v_title


    def getRowsIterator(self, selection = ''):
        """
        Return a row iterator give a certain selection
        Keyword arguments:
        selection -- a selection on the axis of the type "(ant == 'CS001LBA') & (pol == 'XX')"
        """
        if selection != '':
            return self.t.where(selection)
        elif self.selection != '':
            return self.t.where(self.selection)
        else:
            return self.t.iterrows()


    def getValuesGrid(self, selection='', valAxis = "val", notAxes = ["flag"]):
        """
        Try to create a simple matrix of values. NaNs will be returned where the values are not available.
        Keyword arguments:
        selection -- a selection on the axis of the type "(ant == 'CS001LBA') & (pol == 'XX')"
        valAxis -- name of the value axis (use "flag" to obtain the matix of flags)
        notAxis -- list of axes names which are to ignore when looking for all the axes (use "val" when obtaining the matrix of flags) - WARNING: if igoring an axis which index multiple values, then a random value among those possible indexed by that axis is used!
        """

        import numpy as np

        # retreive axes values in a list of tuples
        # [(axisname1,[axisvals1]),(axisname2,[axisvals2]),...]
        axesVals = []
        for axis in self.getAxes(notAxes = notAxes+[valAxis]):
            axesVals.append((axis, np.unique(np.array([x[axis] for x in self.getRowsIterator(selection)]))))
        
        # create an ndarray and fill it with NaNs
        vals = np.ndarray(shape = [len(axis[1]) for axis in axesVals])
        vals[:] = np.NAN
        # refill the array with the correct values when they are available
        for row in self.getRowsIterator(selection):
            pos = []
            for axis in axesVals:
                pos.append(np.where(axis[1]==row[axis[0]])[0])
            vals[pos] = row[valAxis]

        return vals, axesVals


    def getAxes(self, notAxes = ["val","flag"]):
        """
        Return a list with all the axis names in the correct order for
        slicing the getValuesGrid() reurned list.
        Keyword arguments:
        notAxes -- array of names of axes that are not to list
        """
        cols = list(self.t.colpathnames)
        for notAxis in notAxes:
            cols.remove(notAxis)
        return cols
