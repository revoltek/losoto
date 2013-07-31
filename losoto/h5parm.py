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
    def __init__(self, h5parmFile, readonly = True):
        if os.path.isfile(h5parmFile):
            if readonly:
                self.h5parm = tables.openFile(h5parmFile, 'r')
            else:
                self.h5parm = tables.openFile(h5parmFile, 'a')
        else:
            if readonly:
                raise Exception('Missing file '+h5parmFile+'.')
            else:
                # add a compression filter
                f = tables.Filters(complevel=9, complib='zlib')
                self.h5parm = tables.openFile(h5parmFile, filters=f, mode='w')
        
        # if the file is new add the version of the h5parm
        # in losoto._version.__h5parmVersion__


    def makeSolset(self, solsetName = ''):
        """
        Create a new solset, if the provided name is not given or exists
        then it falls back on the first available sol###
        """

        if solsetName in self.getSolsets():
            logging.error('Solution set '+solsetName+' already present.')
            solsetName = ''


        if solsetName == '':
            solsetName = self._fisrtAvailSolsetName()
        
        logging.info('Creating new solution-set '+solsetName+'.')
        return self.h5parm.create_group("/", solsetName)


    def getSolsets(self):
        """
        Return a list of all the available solultion-sets
        """
        return self.h5parm.root._v_children.keys()


    def _fisrtAvailSolsetName(self):
        """
        Create and return the first available solset name which
        has the form of "sol###"
        """
        nums = []
        for solset in self.getSolsets():
#            try:
                if solset[0:3] == 'sol':
                    nums.append(int(solset[3:6]))
#            except:
#                pass

        return "sol%03d" % min(list(set(range(1000)) - set(nums)))


        # masked tables for flags implemented as http://pytables.github.io/cookbook/custom_data_types.html
    
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

        return self.h5parm.createTable(solset, soltabName, descriptor, soltype)


    def getSoltab(self, solset=None):
        """
        Return a list of all the available solultion-tables into a specified solution-set
        """
        if solset == None:
            raise Exception("Solution set not specified while querying for solution-tables list.")

        return solset._v_children.keys()


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
        for soltab in self.getSoltab(solset):
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


def getSols( H, solType, station=[], pol=[], source=[] ):
    """
    Return solution array once selected for solType, station, pol and source.
    solType must be specified.
    """

    if solType == 'amp':
        pass

    elif solType == 'phase':
        pass

    elif solType == 'clock':
        pass

    elif solType == 'TEC':
        pass

    else:
        print "ERROR: unknown solType", solType, "."
        sys.exit(1)

