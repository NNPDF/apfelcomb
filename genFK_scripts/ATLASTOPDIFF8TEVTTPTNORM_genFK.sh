# Merge together the ttbar total xsec grids N_dat times
FKselfmerge FK_ATLASTOPDIFF8TEVTOTDEN.dat 6 > FK_ATLASTOPDIFF8TEVTOTTTPT.dat

# Remove the old, single-datapoint total xsec grid
rm FK_ATLASTOPDIFF8TEVTOTDEN.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_ATLASTOPDIFF8TEVTTPTNUM.dat
FK: FK_ATLASTOPDIFF8TEVTOTTTPT.dat
OP: RATIO" > FK_ATLASTOPDIFF8TEVTTPTNORM-COMPOUND.dat

# Cleanup
exit 0

