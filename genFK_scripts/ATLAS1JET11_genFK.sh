# Run merge
FKmerge FK_ATLAS1JET11_Y1.dat FK_ATLAS1JET11_Y2.dat FK_ATLAS1JET11_Y3.dat FK_ATLAS1JET11_Y4.dat FK_ATLAS1JET11_Y5.dat FK_ATLAS1JET11_Y6.dat > FK_ATLAS1JET11.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_ATLAS1JET11_Y*

exit 0
