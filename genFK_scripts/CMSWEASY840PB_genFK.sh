#!/bin/bash

if [ ! -e "FK_CMSWEASY840PB_WP.dat" ] 
then
	exit 1;
fi

if [ ! -e "FK_CMSWEASY840PB_WM.dat" ]
then
        exit 1;
fi


echo "# COMPOUND FK
FK: FK_CMSWEASY840PB_WP.dat
FK: FK_CMSWEASY840PB_WM.dat
OP: ASY" > FK_CMSWEASY840PB-COMPOUND.dat

# Cleanup

exit 0

