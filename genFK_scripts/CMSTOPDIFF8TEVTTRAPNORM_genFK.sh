# Merge together the ttbar total xsec grids N_dat times
FKselfmerge FK_CMSTOPDIFF8TEVTOTDEN.dat 10 > FK_CMSTOPDIFF8TEVTOTTTRAP.dat

# Remove the old, single-datapoint total xsec grid
rm FK_CMSTOPDIFF8TEVTOTDEN.dat

echo "# COMPOUND FK
FK: FK_CMSTOPDIFF8TEVTTRAPNUM.dat
FK: FK_CMSTOPDIFF8TEVTOTTTRAP.dat
OP: RATIO" > FK_CMSTOPDIFF8TEVTTRAPNORM-COMPOUND.dat

# Cleanup
exit 0

