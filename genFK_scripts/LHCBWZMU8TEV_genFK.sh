# merge FK tables
FKmerge FK_LHCBWZMU8TEV_Z.dat FK_LHCBWZMU8TEV_WP.dat FK_LHCBWZMU8TEV_WM.dat > FK_LHCBWZMU8TEV.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
rm FK_LHCBWZMU8TEV_Z.dat FK_LHCBWZMU8TEV_WP.dat FK_LHCBWZMU8TEV_WM.dat

exit 0
