# Run merge
FKincrement FK_CMSZDIFF12-BIN1.dat 2 > FK_CMSZDIFF12-BIN1.tmp.dat
FKincrement FK_CMSZDIFF12-BIN2.dat 2 > FK_CMSZDIFF12-BIN2.tmp.dat
FKincrement FK_CMSZDIFF12-BIN3.dat 2 > FK_CMSZDIFF12-BIN3.tmp.dat
FKincrement FK_CMSZDIFF12-BIN4.dat 2 > FK_CMSZDIFF12-BIN4.tmp.dat
FKincrement FK_CMSZDIFF12-BIN5.dat 2 > FK_CMSZDIFF12-BIN5.tmp.dat

FKmerge FK_CMSZDIFF12-BIN1.tmp.dat FK_CMSZDIFF12-BIN2.tmp.dat FK_CMSZDIFF12-BIN3.tmp.dat FK_CMSZDIFF12-BIN4.tmp.dat FK_CMSZDIFF12-BIN5.tmp.dat > FK_CMSZDIFF12.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_CMSZDIFF12-BIN*

exit 0
