# merge FK tables
FKmerge FK_CMSWMU8TEV_WP_leptrap.dat FK_CMSWMU8TEV_WM_leptrap.dat > FK_CMSWMU8TEV.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_CMSWMU8TEV_WP_leptrap.dat FK_CMSWMU8TEV_WM_leptrap.dat
exit 0
