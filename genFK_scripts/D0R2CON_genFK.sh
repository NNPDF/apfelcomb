# Run merge
mergeFK FK_D0R2CON_0.dat  FK_D0R2CON_1.dat  FK_D0R2CON_2.dat  FK_D0R2CON_3.dat  FK_D0R2CON_4.dat FK_D0R2CON_5.dat > FK_D0R2CON.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_D0R2CON_*

exit 0
