#!/bin/bash
TH=$1

mkdir err
for i in `seq 1 43`;
do
    FILENAME=discomb$i-$TH
    echo "./dis_comb $i $TH
rc=\$?; if [[ \$rc != 0 ]]; then echo \$rc > ./err/$FILENAME; fi" > $FILENAME
    chmod +x $FILENAME
    echo $FILENAME
    addqueue -m 2 -c "1-12 hrs" $FILENAME
done  
