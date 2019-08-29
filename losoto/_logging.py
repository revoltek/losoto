#!/usr/bin/env python
# encoding: utf-8
import logging

def add_coloring_to_emit_ansi(fn):
    """
    colorize the logging output
    """
    # add methods we need to the class
    def new(*args):
        levelno = args[1].levelno
        if(levelno>=50):
            color = '\x1b[31m' # red
        elif(levelno>=40):
            color = '\x1b[31m' # red
        elif(levelno>=30):
            color = '\x1b[33m' # yellow
        elif(levelno>=20):
            color = '\x1b[32m' # green
        elif(levelno>=10):
            color = '\x1b[35m' # pink
        else:
            color = '\x1b[0m' # normal
        args[1].msg = color + args[1].msg +  '\x1b[0m'  # normal
        return fn(*args)
    return new

# set the logging colors
logging.StreamHandler.emit = add_coloring_to_emit_ansi(logging.StreamHandler.emit)
# set the logging format and default level (info)
logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S", level=logging.INFO)


def setLevel(level):
    """
    Change verbosity
    """
    if level == 'warning':
        logging.root.setLevel(logging.WARNING)
    elif level == 'info':
        logging.root.setLevel(logging.INFO)
    elif level == 'debug':
        logging.root.setLevel(logging.DEBUG)
    else:
        print("Debug level %s doesn't exist." % level)

