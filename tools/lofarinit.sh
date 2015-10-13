#!/bin/sh
# optimized for cep3
if [ -z $PYTHONPATH ]; then export $PYTHONPATH="" ; fi
export PYTHONPATH="/home/fdg/scripts/:/home/fdg/opt/lib/python:/home/fdg/opt/lib/python2.6/site-packages:/usr/lib/python2.6/site-packages:/home/fdg/opt/lib/python2.7/site-packages:/usr/lib/python2.7/site-packages:/home/fdg/.local/lib/python2.7/site-packages:${PYTHONPATH}"
export PATH="${PATH}:/home/fdg/scripts/losoto/bin"
