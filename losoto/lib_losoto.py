#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Some utilities for operations

import os, sys, ast, re
import functools, warnings # for deprecated_alias decorator
from configparser import ConfigParser
if (sys.version_info > (3, 0)):
    #from configparser import ConfigParser
    from io import StringIO
else:
    #from ConfigParser import ConfigParser
    from StringIO import StringIO

import numpy as np
from losoto._logging import logger as logging

cacheSteps = ['plot','clip','flag','norm','smooth'] # steps to use chaced data

class LosotoParser(ConfigParser):
    """
    A parser for losoto parset files.

    Parameters
    ----------
    parsetFile : str
        Name of the parset file.
    """

    def __init__(self, parsetFile):
        ConfigParser.__init__(self, inline_comment_prefixes=('#',';'))

        config = StringIO()
        # add [_global] fake section at beginning
        config.write('[_global]\n'+open(parsetFile).read())
        config.seek(0, os.SEEK_SET)
        self.readfp(config)

    def checkSpelling(self, s, soltab, availValues=[]):
        """
        check if any value in the step is missing from a value list and return a warning
        """
        entries = [x.lower() for x in list(dict(self.items(s)).keys())]
        availValues = ['soltab','operation'] + availValues + \
                    soltab.getAxesNames() + [a+'.minmaxstep' for a in soltab.getAxesNames()] + [a+'.regexp' for a in soltab.getAxesNames()]
        availValues = [x.lower() for x in availValues]
        for e in entries:
            if e not in availValues:
                logging.warning('Mispelled option: %s - Ignoring!' % e)

    def getstr(self, s, v, default=None):
        if self.has_option(s, v):
            return str(self.get(s, v).replace('\'','').replace('"','')) # remove apex
        elif default is None:
            logging.error('Section: %s - Values: %s: required (expected string).' % (s, v))
        else:
            return default

    def getbool(self, s, v, default=None):
        if self.has_option(s, v):
            return ConfigParser.getboolean(self, s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required (expected bool).' % (s, v))
        else:
            return default

    def getfloat(self, s, v, default=None):
        if self.has_option(s, v):
            return ConfigParser.getfloat(self, s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required (expected float).' % (s, v))
        else:
            return default

    def getint(self, s, v, default=None):
        if self.has_option(s, v):
            return ConfigParser.getint(self, s, v)
        elif default is None:
            logging.error('Section: %s - Values: %s: required (expected int).' % (s, v))
        else:
            return default

    def getarray(self, s, v, default=None):
        if self.has_option(s, v):
            try:
                # TODO: why are square brackets being replaced here? What if my selection is a string containing square brackets?
                # return self.getstr(s, v).replace(' ','').replace('[','').replace(']','').split(',') # split also turns str into 1-element lists
                parm = self.getstr(s, v) # split also turns str into 1-element lists
                if parm[0] == '[': # hardcoded for square brackets in parameter set...
                    parm = parm[1:]
                if parm[-1] == ']':
                    parm = parm[:-1]
                return parm.replace(' ','').split(',')
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

    def getarrayfloat2d(self, s, v, default=None):
        """Alternative to parse 1,2 or ndim array input.
           'getarrayfloat() does not support more than 1 dim.
           TODO: it might be cleaner to unify these functions..."""
        try:
            # Remove space after [
            x = self.getstr(s, v, default)
            x = re.sub('\[ +', '[', x.strip())
            # Replace commas and spaces
            x = re.sub('[,\s]+', ', ', x)
            return np.array(ast.literal_eval(x))
        except:
            logging.error('Error interpreting section: %s - values: %s (expected array of float.)' % (s, v))


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
    Return a list of soltabs object for a step and apply selection creteria

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

# fancy backwards compatibility of keywords: allow aliases
# https://stackoverflow.com/questions/49802412/how-to-implement-deprecation-in-python-with-argument-alias#
def deprecated_alias(**aliases):
    def deco(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            rename_kwargs(f.__name__, kwargs, aliases)
            return f(*args, **kwargs)
        return wrapper
    return deco

def rename_kwargs(func_name, kwargs, aliases):
    for alias, new in aliases.items():
        if alias in kwargs:
            if new in kwargs:
                raise TypeError('{} received both {} and {}'.format(
                    func_name, alias, new))
            warnings.warn('{} is deprecated; use {}'.format(alias, new),
                          DeprecationWarning)
            kwargs[new] = kwargs.pop(alias)
