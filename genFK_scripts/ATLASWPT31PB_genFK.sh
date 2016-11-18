# Merge together the WP total xsec grids N_dat times
FKselfmerge FK_ATLASWPT31PB_WP_TOT_1.dat 11 > FK_ATLASWPT31PB_WP_TOT.dat
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Merge together the WM total xsec grids N_dat times
FKselfmerge FK_ATLASWPT31PB_WM_TOT_1.dat 11 > FK_ATLASWPT31PB_WM_TOT.dat
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Remove the old, single-datapoint total xsec grids
rm FK_ATLASWPT31PB_WM_TOT_1.dat FK_ATLASWPT31PB_WP_TOT_1.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_ATLASWPT31PB_WP.dat
FK: FK_ATLASWPT31PB_WM.dat
FK: FK_ATLASWPT31PB_WP_TOT.dat
FK: FK_ATLASWPT31PB_WM_TOT.dat
OP: SMN" > FK_ATLASWPT31PB-COMPOUND.dat

# Cleanup
exit 0
