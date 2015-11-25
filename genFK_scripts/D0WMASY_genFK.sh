#!/bin/bash

if [ ! -e "FK_D0WMASY_WP.dat" ] 
then
	exit 1;
fi

if [ ! -e "FK_D0WMASY_WM.dat" ]
then
        exit 1;
fi


echo "# COMPOUND FK
FK: FK_D0WMASY_WP.dat
FK: FK_D0WMASY_WM.dat
OP: ASY" > FK_D0WMASY-COMPOUND.dat

# Cleanup

exit 0

