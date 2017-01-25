# merge FK tables
FKmerge FK_ATLASWZRAP11_WP.dat FK_ATLASWZRAP11_WM.dat FK_ATLASWZRAP11_Z.dat > FK_ATLASWZRAP11.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_ATLASWZRAP11_WP.dat FK_ATLASWZRAP11_WM.dat FK_ATLASWZRAP11_Z.dat
exit 0
