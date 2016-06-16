#!/bin/bash
set -o xtrace
set -e
set -u
NAME="ATLASWZTOT13TEV81PB"
COMPONENTS="FK_${NAME}_WM_tot.dat FK_${NAME}_WP_tot.dat  FK_${NAME}_Z_tot.dat "
FKmerge $COMPONENTS > FK_${NAME}.dat
sed "s%ATLAS TOTAL W MINUS XS%ATLAS TOTAL W and Z cross sections 13 Tev%g" -i FK_${NAME}.dat
rm ${COMPONENTS}

