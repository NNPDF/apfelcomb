# Run merge
FKmerge FK_CMSZDIFF12-BIN1.dat FK_CMSZDIFF12-BIN2.dat FK_CMSZDIFF12-BIN3.dat FK_CMSZDIFF12-BIN4.dat FK_CMSZDIFF12-BIN5.dat > FK_CMSZDIFF12.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

exit 0
