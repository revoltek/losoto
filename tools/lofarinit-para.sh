#!/bin/sh
# optimized for leiden paracluster
if [ -z $PYTHONPATH ]; then export $PYTHONPATH="" ; fi
export PYTHONPATH="/net/spaarne/data2/scripts:/net/spaarne/data2/scripts/losoto:${PYTHONPATH}"
export PATH="${PATH}:/net/spaarne/data2/scripts/losoto/bin"
