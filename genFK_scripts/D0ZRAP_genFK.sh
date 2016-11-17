#!/bin/bash 

if [ ! -e "FK_D0ZRAP.dat" ] 
then
	exit 1;
fi

FKselfmerge FK_D0ZRAP_TOT_1.dat 28 > FK_D0ZRAP_TOT.dat

RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

echo "# COMPOUND FK
FK: FK_D0ZRAP.dat
FK: FK_D0ZRAP_TOT.dat
OP: RATIO" > FK_D0ZRAP-COMPOUND.dat

# Cleanup
rm FK_D0ZRAP_TOT_1.dat
exit 0
