#!/bin/sh
# optimized for cep3
if [ -z $PYTHONPATH ]; then export $PYTHONPATH="" ; fi
export PYTHONPATH=":/home/fdg/losoto:$PYTHONPATH"
export PATH="/home/fdg/losoto/bin:$PATH"
