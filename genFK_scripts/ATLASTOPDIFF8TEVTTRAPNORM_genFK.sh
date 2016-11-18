# Merge together the ttbar total xsec grids N_dat times
FKselfmerge FK_ATLASTOPDIFF8TEVTOTDEN.dat 10 > FK_ATLASTOPDIFF8TEVTOTTTRAP.dat

# Remove the old, single-datapoint total xsec grid
rm FK_ATLASTOPDIFF8TEVTOTDEN.dat

echo "# COMPOUND FK
FK: FK_ATLASTOPDIFF8TEVTTRAPNUM.dat
FK: FK_ATLASTOPDIFF8TEVTOTTTRAP.dat
OP: RATIO" > FK_ATLASTOPDIFF8TEVTTRAPNORM-COMPOUND.dat

# Cleanup
exit 0

