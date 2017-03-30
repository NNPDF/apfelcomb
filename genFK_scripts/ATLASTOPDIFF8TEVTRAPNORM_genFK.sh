# Merge together the ttbar total xsec grids N_dat times
FKselfmerge FK_ATLASTOPDIFF8TEVTOTDEN.dat 10 > FK_ATLASTOPDIFF8TEVTOTTRAP.dat

# Remove the old, single-datapoint total xsec grid
rm FK_ATLASTOPDIFF8TEVTOTDEN.dat

# write the compound file
echo "# COMPOUND FK
FK: FK_ATLASTOPDIFF8TEVTRAPNUM.dat
FK: FK_ATLASTOPDIFF8TEVTOTTRAP.dat
OP: RATIO" > FK_ATLASTOPDIFF8TEVTRAPNORM-COMPOUND.dat

# Cleanup
exit 0
