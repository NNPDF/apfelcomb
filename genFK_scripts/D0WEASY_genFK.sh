#!/bin/bash

if [ ! -e "FK_D0WEASY_WP.dat" ] 
then
	exit 1;
fi

if [ ! -e "FK_D0WEASY_WM.dat" ]
then
        exit 1;
fi


echo "# COMPOUND FK
FK: FK_D0WEASY_WP.dat
FK: FK_D0WEASY_WM.dat
OP: ASY" > FK_D0WEASY-COMPOUND.dat

# Cleanup

exit 0

