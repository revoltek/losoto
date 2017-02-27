#!/bin/sh
# optimized for leiden paracluster
if [ -z $PYTHONPATH ]; then export $PYTHONPATH="" ; fi
export PYTHONPATH="/net/lofar1/data1/oonk/rh7_losoto_generic:${PYTHONPATH}"
export PATH="${PATH}:/net/lofar1/data1/oonk/rh7_losoto_generic/bin"
