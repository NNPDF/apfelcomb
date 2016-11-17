# Merge together the ttbar total xsec grids N_dat times
FKselfmerge FK_ATLASTOPDIFF8TEVTOTDEN.dat 7 > FK_ATLASTOPDIFF8TEVTOTTTM.dat

# Remove the old, single-datapoint total xsec grid
rm FK_ATLASTOPDIFF8TEVTOTDEN.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_ATLASTOPDIFF8TEVTTMNUM.dat
FK: FK_ATLASTOPDIFF8TEVTOTTTM.dat
OP: RATIO" > FK_ATLASTOPDIFF8TEVTTMNORM-COMPOUND.dat

# Cleanup

exit 0

