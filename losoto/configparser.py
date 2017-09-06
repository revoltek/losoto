"""
ConfigFile parser
"""

import logging
import ConfigParser
import json

def parset(filename):
    """
    Parset parser
    """
    logging.debug("Parsing %s" % filename)
    parser = ConfigParser.MagicConfigParser()
    data = parser.read(filename)
    # create a dict for each keyword
    p = {}
    # iterate on operations
    for step in :
        p[step] = {}
        # iterate over keyword
        for keyword in :
            p[step][keyword] = json.loads(config.get(step,keyword))

    return p
