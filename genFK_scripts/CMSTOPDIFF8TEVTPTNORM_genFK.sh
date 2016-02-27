#!/bin/bash

if [ ! -e "FK_CMSTOPDIFF8TEVTPT.dat" ] 
then
	exit 1;
fi

if [ ! -e "FK_CMSTOPDIFF8TEVTOT.dat" ]
then
        exit 1;
fi


echo "# COMPOUND FK
FK: FK_CMSTOPDIFF8TEVTPT.dat
FK: FK_CMSTOPDIFF8TEVTOT.dat
OP: RATIO" > FK_CMSTOPDIFF8TEVTPTNORM-COMPOUND.dat

# Cleanup

exit 0

