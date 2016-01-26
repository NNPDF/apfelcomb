#!/bin/bash
TH=$1

mkdir err
for i in `seq 1 37`;
do
    FILENAME=discomb$i-$TH
    echo "./dis_comb $i $TH
rc=\$?; if [[ \$rc != 0 ]]; then echo \$rc > ./err/$FILENAME; fi" > $FILENAME
    chmod +x $FILENAME
    echo $FILENAME
    mpisubnopause "1-9 hrs" 1 $FILENAME
done  
