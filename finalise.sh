#!/bin/bash
# This script finalises a theory

echo "Finalising theoryID: " $1
#./run_cfactors.py $1
cp -r ./compound ./results/theory_$1/
rm -rf ./results/theory_$1/subgrids/
cd results
tar -cvzf ./theory_$1.tgz ./theory_$1/
cd ../
