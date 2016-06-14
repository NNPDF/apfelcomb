# merge FK tables
FKmerge FK_CMSTTBARTOT7TEV.dat CMSTTBARTOT8TEV.dat > FK_CMSTTBARTOT.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
rm FK_CMSTTBARTOT7TEV.dat CMSTTBARTOT8TEV.dat

exit 0
