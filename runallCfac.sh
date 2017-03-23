#!/bin/bash

dataDir="../nnpdfcpp/data/"
sourceDir=$dataDir"NNLOCFAC/"
cFactors=$sourceDir"registeredCFactors.dat"

echo "Processing theoryID: " $1

while read p; do
  a=( $p )
  if [ "${a[0]}" = "QCD" ]; then
  	echo ${a[0]} ${a[1]} ${a[2]} ${a[3]}
  	./src/cfac_scale $1 ${a[1]} ${a[2]} ${a[3]}
  else
  	echo cp $sourceDir"CF_"${a[0]}"_"${a[1]}".dat" $dataDir"theory_"$1"/cfactor/"
  	cp $sourceDir"CF_"${a[0]}"_"${a[1]}".dat" $dataDir"theory_"$1"/cfactor/"
  fi
done < $cFactors