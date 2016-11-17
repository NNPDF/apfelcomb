# Merge together the ttbar total xsec grids N_dat times
FKselfmerge FK_CMSTOPDIFF8TEVTOTDEN.dat 8 > FK_CMSTOPDIFF8TEVTOTTPT.dat

# Remove the old, single-datapoint total xsec grid
rm FK_CMSTOPDIFF8TEVTOTDEN.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_CMSTOPDIFF8TEVTPTNUM.dat
FK: FK_CMSTOPDIFF8TEVTOTTPT.dat
OP: RATIO" > FK_CMSTOPDIFF8TEVTPTNORM-COMPOUND.dat

# Cleanup
exit 0

