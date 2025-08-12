#!/bin/bash
#
# Compare the outputs generated with benchmark_all.sh

[ $# -eq 2 ] || { echo "Usage <reference-dir> <output-dir>"; exit 1; }

ref_dir=$1
out_dir=$2

for operation in flag flagextend flagstation plot reweight
do
  echo -e "\n\n**** Comparing reference and output for operation ${operation^^} ..."
  echo -e "\nDifferences between log files:"
  diff --side-by-side --suppress-common-lines \
    <(grep -v RuntimeWarning ${ref_dir}/${operation}.log | cut -b23- | sort) \
    <(grep -v RuntimeWarning ${out_dir}/${operation}.log | cut -b23- | sort)
  if [ "${operation}" = "plot" ]
  then
    echo -e "\nDifferences between generated plots:"
    for ref in ${ref_dir}/*.png
    do
      out=${out_dir}/${ref##*/}
      diff -q ${ref} ${out}
    done
  else
    echo -e "\nDifferences between H5 modifications:"
    diff --side-by-side --suppress-common-lines \
      ${ref_dir}/h5diff.${operation}.txt \
      ${out_dir}/h5diff.${operation}.txt
  fi
  echo -e "\nDifferences in execution times:"
  diff --side-by-side --suppress-common-lines --width 60 \
    ${ref_dir}/${operation}.time \
    ${out_dir}/${operation}.time
done
