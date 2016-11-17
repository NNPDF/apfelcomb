# Merge together the ttbar total xsec grids N_dat times
FKselfmerge FK_ATLASTOPDIFF8TEVTOTDEN.dat 8 > FK_ATLASTOPDIFF8TEVTOTTPT.dat

# Remove the old, single-datapoint total xsec grid
rm FK_ATLASTOPDIFF8TEVTOTDEN.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_ATLASTOPDIFF8TEVTPTNUM.dat
FK: FK_ATLASTOPDIFF8TEVTOTTPT.dat
OP: RATIO" > FK_ATLASTOPDIFF8TEVTPTNORM-COMPOUND.dat

# Cleanup
exit 0

