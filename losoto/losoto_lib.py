#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Some utilities for operations

import os, sys
import ast, re
import logging

class LosotoParser(object):

    def __init__(self, parsetFile):
        import StringIO, ConfigParser
        self.parsetFile = parsetFile
        config = StringIO.StringIO()
        config.write('[_global]\n'+open(parsetFile).read().replace('#',';'))
        config.seek(0, os.SEEK_SET)
        self.parser = ConfigParser.RawConfigParser()
        self.parser.readfp(config)

    def sections(self):
        return self.parser.sections()

    def getstr(self, s, v, default=None):
        if self.parser.has_option(s, v):
            return self.parser.get(s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required.' % (s, v))
        else:
            return default

    def getbool(self, s, v, default=None):
        if self.parser.has_option(s, v):
            return self.parser.getboolean(s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required.' % (s, v))
        else:
            return default

    def getfloat(self, s, v, default=None):
        if self.parser.has_option(s, v):
            return self.parser.getfloat(s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required.' % (s, v))
        else:
            return default

    def getint(self, s, v, default=None):
        if self.parser.has_option(s, v):
            return self.parser.getint(s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required.' % (s, v))
        else:
            return default

    def getarray(self, s, v, default=None):
        if self.parser.has_option(s, v):
            try:
                return ast.literal_eval( self.parser.get(s, v) )
            except:
                logging.error('Problem interpreting section: %s - values: %s' % (s, v))
                sys.exit(1)
        elif default is None:
            logging.error('Section: %s - Values: %s: required.' % (s, v))
        else:
            return default

    def has_option(self, s, v):
        return self.parser.has_option(s, v)


def getParAxis( parser, step, axisName ):
    """
    Parameters
    ----------
    parser : parser obj
        configuration file
    step : str
        this step
    axisName : str
        an axis name

    Returns
    -------
    str, dict or list
        a selection criteria
    """
    axisOpt = None
    if parser.has_option(step, axisName):
        axisOpt = parser.getstr(step, axisName)
        # if vector/dict, reread it
        if axisOpt != '' and (axisOpt[0] == '[' and axisOpt[-1] == ']') or (axisOpt[0] == '{' and axisOpt[-1] == '}'):
            axisOpt = ast.literal_eval( parser.getstr(step, axisName) )

    elif parser.has_option('_global', axisName):
        axisOpt = parser.getstr('_global', axisName)
        # if vector/dict, reread it
        if axisOpt != '' and (axisOpt[0] == '[' and axisOpt[-1] == ']') or (axisOpt[0] == '{' and axisOpt[-1] == '}'):
            axisOpt = ast.literal_eval( parser.getstr('_global', axisName) )

    if axisOpt == '' or axisOpt == []:
        axisOpt = None
 
    return axisOpt


def getStepSoltabs(parser, step, H):
    """
    Return a list of soltabs object for 

    Parameters
    ----------
    parser : parser obj
        configuration file
    step : str
        current step
    H : h5parm obj
        the h5parm object

    Returns
    -------
    list
        list of soltab obj with applied selection
    """
    cacheSteps = ['clip','flag'] # steps to use chaced data

    # selection on soltabs
    if parser.has_option(step, 'soltab'):
        stsel = ast.literal_eval(parser.getstr(step, 'soltab'))
    elif parser.has_option('_global', 'soltab'):
        stsel = ast.literal_eval(parser.getstr('_global', 'soltab'))
    else:
        stsel = '*/*' # select all
    if not type(stsel) is list: stsel = [stsel]

    soltabs = []
    for solset in H.getSolsets():
        for soltabName in solset.getSoltabNames():
            if any(re.compile(this_stsel).match(solset.name+'/'+soltabName) for this_stsel in stsel):
                if step in cacheSteps:
                    soltabs.append( solset.getSoltab(soltabName, useCache=True) )
                else:
                    soltabs.append( solset.getSoltab(soltabName, useCache=False) )

    # axes selection
    for soltab in soltabs:
        userSel = {}
        for axisName in soltab.getAxesNames():
            userSel[axisName] = getParAxis( parser, step, axisName )
        soltab.setSelection(**userSel)

    return soltabs
