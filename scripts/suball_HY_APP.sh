#!/bin/bash
TH=$1

for i in `seq 0 76`;
do
    FILENAME=apfelcomb$i-$TH
    echo "../appl_comb $i $TH" > $FILENAME
    chmod +x $FILENAME
    echo $FILENAME
    mpisubnopause "1-9 hrs" 1 $FILENAME
done  
