#!/bin/bash
if [ ! -e "FK_CMSWCHARM_WP.dat" ]
then
	exit 1;
fi

if [ ! -e "FK_CMSWCHARM_WM.dat" ]
then
	exit 1;
fi

echo "# COMPOUND FK
FK: FK_CMSWCHARM_WP.dat
FK: FK_CMSWCHARM_WM.dat
OP: RATIO" > FK_CMSWCHARMRAT-COMPOUND.dat

echo "# COMPOUND FK
FK: FK_CMSWCHARM_WP.dat
FK: FK_CMSWCHARM_WM.dat
OP: ADD" > FK_CMSWCHARMTOT-COMPOUND.dat

# Cleanup

exit 0
