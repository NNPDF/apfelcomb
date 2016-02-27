#!/bin/bash

if [ ! -e "FK_CMSTOPDIFF8TEVTTPT.dat" ] 
then
	exit 1;
fi

if [ ! -e "FK_CMSTOPDIFF8TEVTOT.dat" ]
then
        exit 1;
fi


echo "# COMPOUND FK
FK: FK_CMSTOPDIFF8TEVTTPT.dat
FK: FK_CMSTOPDIFF8TEVTOT.dat
OP: RATIO" > FK_CMSTOPDIFF8TEVTTPTNORM-COMPOUND.dat

# Cleanup

exit 0

