#!/bin/bash
TH=$1

for i in `seq 0 143`;
do
    FILENAME=apfelcomb$i-$TH
    echo "./appl_comb $i $TH
rc=\$?; if [[ \$rc != 0 ]]; then echo \$rc > ./err/$FILENAME; fi" > $FILENAME
    chmod +x $FILENAME
    echo $FILENAME
    addqueue -c "1-12 hrs" $FILENAME
done  
