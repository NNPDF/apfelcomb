#!/bin/bash
# This script finalises a theory

echo "Processing theoryID: " $1
./merge_allgrids.py $1
./run_cfactors.py $1
cp -r ./compound ./results/theory_$1/