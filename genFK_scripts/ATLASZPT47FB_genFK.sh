# Merge together the WP total xsec grids N_dat=26 times
FKselfmerge FK_ATLASZPT47FB_TOT_1.dat 26 >  FK_ATLASZPT47FB_TOT.dat
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Remove the old, single-datapoint total xsec grid
rm FK_ATLASZPT47FB_TOT_1.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_ATLASZPT47FB.dat
FK: FK_ATLASZPT47FB_TOT.dat
OP: RATIO" > FK_ATLASZPT47FB-COMPOUND.dat

# Cleanup
exit 0
