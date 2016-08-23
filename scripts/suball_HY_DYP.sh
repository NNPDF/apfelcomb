#!/bin/bash
TH=$1

mkdir err
for i in `seq 1 23`;
do
    FILENAME=dypcomb$i-$TH
    echo "./ftdy_comb $i $TH
rc=\$?; if [[ \$rc != 0 ]]; then echo \$rc > ./err/$FILENAME; fi" > $FILENAME
    chmod +x $FILENAME
    echo $FILENAME
    addqueue -c "1-12 hrs" $FILENAME
done  
