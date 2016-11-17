# Merge together the ttbar total xsec grids N_dat times
FKmerge FK_CMSTOPDIFF8TEVTOTDEN.dat 6 > FK_CMSTOPDIFF8TEVTOTTTPT.dat

# Remove the old, single-datapoint total xsec grid
rm FK_CMSTOPDIFF8TEVTOTDEN.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_CMSTOPDIFF8TEVTTPTNUM.dat
FK: FK_CMSTOPDIFF8TEVTOTTTPT.dat
OP: RATIO" > FK_CMSTOPDIFF8TEVTTPTNORM-COMPOUND.dat

# Cleanup
exit 0

