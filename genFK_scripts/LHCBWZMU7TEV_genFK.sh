# merge FK tables
FKmerge FK_LHCBWZMU7TEV_Z.dat FK_LHCBWZMU7TEV_WP.dat FK_LHCBWZMU7TEV_WM.dat > FK_LHCBWZMU7TEV.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
rm FK_LHCBWZMU7TEV_Z.dat FK_LHCBWZMU7TEV_WP.dat FK_LHCBWZMU7TEV_WM.dat

exit 0
