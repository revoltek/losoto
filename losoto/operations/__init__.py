import os, time, glob
from os.path import dirname, basename, isfile

modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py') and not basename(f).startswith('_')]
del modules
from . import *

class Timer(object):
    """
    context manager used to time the operations
    """

    def __init__(self, log=None, step = 'undef.', operation = 'undef.'):
        """
        log: is a logging istance to print the correct log format
        if nothing is passed, root is used
        """
        if log is None: self.log = logging
        else: self.log = log
        self.step = step
        self.operation = operation

    def __enter__(self):
        self.log.info("--> Starting \'" + self.step + "\' step (operation: " + self.operation + ").")
        self.start = time.time()
        self.startcpu = time.perf_counter()

    def __exit__(self, exit_type, value, tb):

        # if not an error
        if exit_type is None:
            self.log.info("Time for %s step: %i s (cpu: %i s)." % ( self.step, ( time.time() - self.start), (time.perf_counter() - self.startcpu) ))
