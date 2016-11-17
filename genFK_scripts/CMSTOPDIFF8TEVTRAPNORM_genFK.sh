# Merge together the ttbar total xsec grids N_dat times
FKselfmerge FK_CMSTOPDIFF8TEVTOTDEN.dat 10 > FK_CMSTOPDIFF8TEVTOTTRAP.dat

# Remove the old, single-datapoint total xsec grid
rm FK_CMSTOPDIFF8TEVTOTDEN.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_CMSTOPDIFF8TEVTRAPNUM.dat
FK: FK_CMSTOPDIFF8TEVTOTTRAP.dat
OP: RATIO" > FK_CMSTOPDIFF8TEVTRAPNORM-COMPOUND.dat

# Cleanup

exit 0

