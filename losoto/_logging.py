import os, sys, logging, time

class _ColorStreamHandler(logging.StreamHandler):

    DEFAULT = '\x1b[0m'
    RED     = '\x1b[31m'
    GREEN   = '\x1b[32m'
    YELLOW  = '\x1b[33m'
    PURP    = '\x1b[35m'

    CRITICAL = RED
    ERROR    = RED
    WARNING  = YELLOW
    INFO     = GREEN
    DEBUG    = PURP

    @classmethod
    def _get_color(cls, level):
        if level >= logging.CRITICAL:  return cls.CRITICAL
        elif level >= logging.ERROR:   return cls.ERROR
        elif level >= logging.WARNING: return cls.WARNING
        elif level >= logging.INFO:    return cls.INFO
        elif level >= logging.DEBUG:   return cls.DEBUG
        else:                          return cls.DEFAULT

    def __init__(self, stream=None):
        logging.StreamHandler.__init__(self, stream)

    def format(self, record):
        color = self._get_color(record.levelno)
        record.msg = color + record.msg + self.DEFAULT
        return logging.StreamHandler.format(self, record)

class Logger():

    def __init__(self, level = 'info', logfile = None, log_dir = None):

        self.logfile = logfile
        self.log_dir = log_dir
        self.backup(logfile, log_dir)
        self.set_logger(logfile, log_dir)
        self.set_level(level)

    def backup(self, logfile, log_dir):

        if self.logfile is not None:
            # bkp old log dir
            if os.path.isdir(log_dir):
                current_time = time.localtime()
                log_dir_old = time.strftime(log_dir+'_bkp_%Y-%m-%d_%H:%M', current_time)
                os.system('mv %s %s' % ( log_dir, log_dir_old ))
            os.makedirs(log_dir)
    
            # bkp old log file
            if os.path.exists(logfile):
                current_time = time.localtime()
                logfile_old = time.strftime(logfile+'_bkp_%Y-%m-%d_%H:%M', current_time)
                os.system('mv %s %s' % ( logfile, logfile_old ))
            

    def set_logger(self, logfile, log_dir):
      
        self.logger = logging.getLogger("LoSoTo")
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")

        # create file handler which logs even debug messages
        if self.logfile is not None:
            handlerFile = logging.FileHandler(logfile)
            #handlerFile.setLevel(logging.DEBUG)
            handlerFile.setFormatter(formatter)
            self.logger.addHandler(handlerFile)
        
        # create console handler with a higher log level
        handlerConsole = _ColorStreamHandler(stream=sys.stdout)
        #handlerConsole.setLevel(logging.INFO)
        handlerConsole.setFormatter(formatter)
        self.logger.addHandler(handlerConsole)

    def set_level(self, level):
        if level == 'warning':
            self.logger.setLevel(logging.WARNING)
        elif level == 'info':
            self.logger.setLevel(logging.INFO)
        elif level == 'debug':
            self.logger.setLevel(logging.DEBUG)
        else:
            print("Debug level %s doesn't exist." % level)

# this is used by all libraries for logging
logger = logging.getLogger("LoSoTo")
