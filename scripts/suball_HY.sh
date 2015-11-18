#!/bin/bash
TH=21

for i in `seq 0 70`;
do
    FILENAME=apfelcomb$i-$TH
    echo "./appl_comb $i $TH
./dis_comb $i $TH
./ftdy_comb $i $TH" > $FILENAME
    chmod +x $FILENAME
    echo $FILENAME
    mpisubnopause "1-9 hrs" 1 $FILENAME
done  
