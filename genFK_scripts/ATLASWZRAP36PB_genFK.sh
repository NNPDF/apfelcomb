# merge FK tables
FKmerge FK_ATLAS_WP_2010_36pb.dat FK_ATLAS_WM_2010_36pb.dat FK_ATLAS_Z_2010_36pb.dat > FK_ATLASWZRAP36PB.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_ATLAS_WP_2010_36pb.dat FK_ATLAS_WM_2010_36pb.dat FK_ATLAS_Z_2010_36pb.dat

exit 0
