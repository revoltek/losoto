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
        complevel -- compression level from 0 to 9 (default=9)
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
        """
        if solset == None:
            raise Exception("Solution set not specified while adding a solution-table.")
        if soltype == None:
            raise Exception("Solution type not specified while adding a solution-table.")
        
        soltabName = self._fisrtAvailSoltabName(solset, soltype)
        logging.info('Creating new solution-table '+soltabName+'.')

        return self.H.createTable(solset, soltabName, descriptor, soltype)


    def getSoltabs(self, solset=None):
        """
        Return a dict of all the available solultion-tables into a specified solution-set
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
        add a single row to the given soltab
        """
        if soltab == None:
            raise Exception("Solution-table not specified while adding a new row.")

        soltab.append(val)

     def getAnt(self):
         pass

     def getSou(self):
         pass

class solFetcher():

    def __init__(self, table):
        """
        Keyword arguments:
        tab -- table object
        """
        
        self.t = table







