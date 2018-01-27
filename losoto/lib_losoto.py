#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Some utilities for operations

import os, sys, ast, re
import logging
if (sys.version_info > (3, 0)):
    from configparser import RawConfigParser
else:
    from ConfigParser import RawConfigParser

cacheSteps = ['plot','clip','flag','norm','smooth'] # steps to use chaced data

class LosotoParser(RawConfigParser):
    """
    A parser for losoto parset files.

    Parameters
    ----------
    parsetFile : str
        Name of the parset file.
    """

    def __init__(self, parsetFile):
        RawConfigParser.__init__(self)

        # read parset and replace '#' with ';' to allow # as inline comments
        # also add [_global] fake section at beginning
        import StringIO
        config = StringIO.StringIO()
        config.write('[_global]\n'+open(parsetFile).read().replace('#',';'))
        config.seek(0, os.SEEK_SET)
        self.readfp(config)

    def getstr(self, s, v, default=None):
        if self.has_option(s, v):
            return self.get(s, v).replace('\'','').replace('"','') # remove apex
        elif default is None:
            logging.error('Section: %s - Values: %s: required (expected string).' % (s, v))
        else:
            return default

    def getbool(self, s, v, default=None):
        if self.has_option(s, v):
            return RawConfigParser.getboolean(self, s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required (expected bool).' % (s, v))
        else:
            return default

    def getfloat(self, s, v, default=None):
        if self.has_option(s, v):
            return RawConfigParser.getfloat(self, s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required (expected float).' % (s, v))
        else:
            return default

    def getint(self, s, v, default=None):
        if self.has_option(s, v):
            return RawConfigParser.getint(self, s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required (expected int).' % (s, v))
        else:
            return default

    def getarray(self, s, v, default=None):
        if self.has_option(s, v):
            try:
                return self.getstr(s, v).replace(' ','').replace('[','').replace(']','').split(',') # split also turns str into 1-element lists
            except:
                logging.error('Error interpreting section: %s - values: %s (should be a list as [xxx,yyy,zzz...])' % (s, v))
        elif default is None:
            logging.error('Section: %s - Values: %s: required.' % (s, v))
        else:
            return default

    def getarraystr(self, s, v, default=None):
        try:
            return [str(x) for x in self.getarray(s, v, default)]
        except:
            logging.error('Error interpreting section: %s - values: %s (expected array of str.)' % (s, v))

    def getarraybool(self, s, v, default=None):
        try:
            return [bool(x) for x in self.getarray(s, v, default)]
        except:
            logging.error('Error interpreting section: %s - values: %s (expected array of bool.)' % (s, v))

    def getarrayfloat(self, s, v, default=None):
        try:
            return [float(x) for x in self.getarray(s, v, default)]
        except:
            logging.error('Error interpreting section: %s - values: %s (expected array of float.)' % (s, v))

    def getarrayint(self, s, v, default=None):
        try:
            return [int(x) for x in self.getarray(s, v, default)]
        except:
            logging.error('Error interpreting section: %s - values: %s (expected array of int.)' % (s, v))


def getParAxis( parser, step, axisName ):
    """
    Axes values can be selected with:
    - an array of values
    - a string (reg exp): passed as axisName.regexp = ...
    - a min,max,step format: passed as axisName.minmaxstep = ...

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
        axisOpt = parser.getarray(step, axisName)
    elif parser.has_option(step, axisName+'.regexp'):
        axisOpt = parser.getstr(step, axisName+'.regexp')
    elif parser.has_option(step, axisName+'.minmaxstep'):
        axisOpt = parser.getarray(step, axisName+'.minmaxstep')
        if len(axisOpt) == 3:
            axisOpt = {'min':float(axisOpt[0]), 'max':float(axisOpt[1]), 'step':int(axisOpt[2])}
        else:
            axisOpt = {'min':float(axisOpt[0]), 'max':float(axisOpt[1])}

    # global options
    elif parser.has_option('_global', axisName):
        axisOpt = parser.getarray('_global', axisName)
    elif parser.has_option('_global', axisName+'.regexp'):
        axisOpt = parser.getstr('_global', axisName+'.regexp')
    elif parser.has_option('_global', axisName+'.minmaxstep'):
        axisOpt = parser.getarray('_global', axisName+'.minmaxstep')
        if len(axisOpt) == 3:
            axisOpt = {'min':float(axisOpt[0]), 'max':float(axisOpt[1]), 'step':int(axisOpt[2])}
        else:
            axisOpt = {'min':float(axisOpt[0]), 'max':float(axisOpt[1])}

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

    # selection on soltabs
    if parser.has_option(step, 'soltab'):
        stsel = parser.getarraystr(step, 'soltab')
    elif parser.has_option('_global', 'soltab'):
        stsel = parser.getarraystr('_global', 'soltab')
    else:
        stsel = ['.*/.*'] # select all
    #if not type(stsel) is list: stsel = [stsel]

    soltabs = []
    for solset in H.getSolsets():
        for soltabName in solset.getSoltabNames():
            if any(re.compile(this_stsel).match(solset.name+'/'+soltabName) for this_stsel in stsel):
                if parser.getstr(step, 'operation').lower() in cacheSteps:
                    soltabs.append( solset.getSoltab(soltabName, useCache=True) )
                else:
                    soltabs.append( solset.getSoltab(soltabName, useCache=False) )

    if soltabs == []:
        logging.warning('No soltabs selected for step %s.' % step)

    # axes selection
    for soltab in soltabs:
        userSel = {}
        for axisName in soltab.getAxesNames():
            userSel[axisName] = getParAxis( parser, step, axisName )
        soltab.setSelection(**userSel)

    return soltabs
