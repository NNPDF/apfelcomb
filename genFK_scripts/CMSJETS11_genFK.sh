# Run merge
FKmerge FK_CMSJETS11_ETA1.dat FK_CMSJETS11_ETA2.dat FK_CMSJETS11_ETA3.dat FK_CMSJETS11_ETA4.dat FK_CMSJETS11_ETA5.dat > FK_CMSJETS11.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_CMSJETS11_ETA*
exit 0
