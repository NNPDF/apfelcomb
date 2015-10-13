#!/bin/bash

if [ ! -e "FK_CMSWMASY47FB_WP.dat" ] 
then
	exit 1;
fi

if [ ! -e "FK_CMSWMASY47FB_WM.dat" ]
then
        exit 1;
fi

echo "# COMPOUND FK
FK: FK_CMSWMASY47FB_WP.dat
FK: FK_CMSWMASY47FB_WM.dat
OP: ASY" > FK_CMSWMASY47FB-COMPOUND.dat

# Cleanup

exit 0
