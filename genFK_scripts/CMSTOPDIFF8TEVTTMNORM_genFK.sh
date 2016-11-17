# Merge together the ttbar total xsec grids N_dat times
FKmerge FK_CMSTOPDIFF8TEVTOTDEN.dat 7 > FK_CMSTOPDIFF8TEVTOTTTM.dat

# Remove the old, single-datapoint total xsec grid
rm FK_CMSTOPDIFF8TEVTOTDEN.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_CMSTOPDIFF8TEVTTMNUM.dat
FK: FK_CMSTOPDIFF8TEVTOTTTM.dat
OP: RATIO" > FK_CMSTOPDIFF8TEVTTMNORM-COMPOUND.dat

# Cleanup

exit 0

