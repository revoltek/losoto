#!/bin/bash
#
# Benchmark all the LoSoTo operations that use multiprocessing

for operation in flag flagextend flagstation plot reweight
do
  echo -e "\n\n**** Benchmarking operation ${operation^^} ... ****\n"
  cp bandpass.input.h5 bandpass.h5
  (time losoto --verbose bandpass.h5 ${operation}.parset) |& tee ${operation}.log
  h5diff -v bandpass.input.h5 bandpass.h5 > h5diff.${operation}.txt
  tail -3 ${operation}.log > ${operation}.time
done
