#!/bin/bash

dataDir="../nnpdfcpp/data"
cFactors=$dataDir"/NNLOCFAC/registeredCFactors.dat"

echo "Processing theoryID: " $1

while read p; do
  a=( $p )
  echo ${a[0]} ${a[1]}
  ./src/cfac_scale $1 ${a[0]} ${a[1]}
done < $cFactors