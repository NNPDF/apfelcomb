# Run merge
mergeFK FK_CDFR2KT_0.dat  FK_CDFR2KT_1.dat  FK_CDFR2KT_2.dat  FK_CDFR2KT_3.dat  FK_CDFR2KT_4.dat > FK_CDFR2KT.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_CDFR2KT_*

exit 0
