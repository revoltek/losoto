import os, time, glob

__all__ = [ os.path.basename(f)[:-3] for f in glob.glob(os.path.dirname(__file__)+"/*.py") if not f.endswith('__init__.py')]
for x in __all__: __import__(x, locals(), globals())

class timer(object):
    """
    context manager used to time the operations
    """

    def __init__(self, log='', step = None, operation = None):
        """
        log: is a logging istance to print the correct log format
        if nothing is passed, root is used
        """
        if log == '': self.log = logging
        else: self.log = log
        self.step = step
        self.operation = operation

    def __enter__(self):
        self.log.info("--> Starting \'" + self.step + "\' step (operation: " + self.operation + ").")
        self.start = time.time()
        self.startcpu = time.clock()

    def __exit__(self, type, value, tb):

        if type is not None:
            raise type, value, tb
        self.log.info("Time for this step: %i s (cpu: %i s)." % ( ( time.time() - self.start), (time.clock() - self.startcpu) ))

