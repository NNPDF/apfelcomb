# Run merge
FKmerge FK_CMSZDIFF12-BIN1_ptZ.dat FK_CMSZDIFF12-BIN2_ptZ.dat FK_CMSZDIFF12-BIN3_ptZ.dat FK_CMSZDIFF12-BIN4_ptZ.dat FK_CMSZDIFF12-BIN5_ptZ.dat > FK_CMSZDIFF12.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

rm FK_CMSZDIFF12-BIN*
exit 0
