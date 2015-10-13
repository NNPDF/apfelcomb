# Run merge
mergeFK FK_ATLASR04JETS2P76TEV_ETA1.dat FK_ATLASR04JETS2P76TEV_ETA2.dat FK_ATLASR04JETS2P76TEV_ETA3.dat FK_ATLASR04JETS2P76TEV_ETA4.dat FK_ATLASR04JETS2P76TEV_ETA5.dat FK_ATLASR04JETS2P76TEV_ETA6.dat FK_ATLASR04JETS2P76TEV_ETA7.dat > FK_ATLASR04JETS2P76TEV.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_ATLASR04JETS2P76TEV_ETA*

exit 0
