############## bin 1 ################################
## Merge together the Z total xsec grids N_applgrid times
FKmerge FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat  FK_ATLASZPT7TEV-BIN1_tot.dat FK_ATLASZPT7TEV-BIN1_tot.dat >  FK_ATLASZPT7TEV-BIN1_TOT.dat

# Add 12 dummy points below the Z pT cut
FKincrement FK_ATLASZPT7TEV-BIN1_TOT.dat 12 > FKincrement FK_ATLASZPT7TEV-BIN1_TOT.tmp.dat
FKincrement FK_ATLASZPT7TEV-BIN1_ptZ.dat 12 > FKincrement FK_ATLASZPT7TEV-BIN1_ptZ.tmp.dat

############## bin 2 ################################
FKmerge FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat  FK_ATLASZPT7TEV-BIN2_tot.dat FK_ATLASZPT7TEV-BIN2_tot.dat >  FK_ATLASZPT7TEV-BIN2_TOT.dat

FKincrement FK_ATLASZPT7TEV-BIN2_TOT.dat 12 > FKincrement FK_ATLASZPT7TEV-BIN2_TOT.tmp.dat
FKincrement FK_ATLASZPT7TEV-BIN2_ptZ.dat 12 > FKincrement FK_ATLASZPT7TEV-BIN2_ptZ.tmp.dat

############## bin 3 ################################
FKmerge FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat  FK_ATLASZPT7TEV-BIN3_tot.dat FK_ATLASZPT7TEV-BIN3_tot.dat >  FK_ATLASZPT7TEV-BIN3_TOT.dat

FKincrement FK_ATLASZPT7TEV-BIN3_TOT.dat 12 > FKincrement FK_ATLASZPT7TEV-BIN3_TOT.tmp.dat
FKincrement FK_ATLASZPT7TEV-BIN3_ptZ.dat 12 > FKincrement FK_ATLASZPT7TEV-BIN3_ptZ.tmp.dat

# Remove the old, single-datapoint total xsec grid
rm FK_ATLASZPT7TEV-BIN1_tot.dat
rm FK_ATLASZPT7TEV-BIN2_tot.dat
rm FK_ATLASZPT7TEV-BIN3_tot.dat

## Merge the total grids
FKmerge FK_ATLASZPT7TEV-BIN1_TOT.tmp.dat  FK_ATLASZPT7TEV-BIN2_TOT.tmp.dat  FK_ATLASZPT7TEV-BIN3_TOT.tmp.dat >  FK_ATLASZPT7TEV_TOT.dat
FKmerge FK_ATLASZPT7TEV-BIN1_ptZ.tmp.dat  FK_ATLASZPT7TEV-BIN2_ptZ.tmp.dat  FK_ATLASZPT7TEV-BIN3_ptZ.tmp.dat >  FK_ATLASZPT7TEV_ptZ.dat

## write the compound file for each bin
echo "# COMPOUND FK
FK: FK_ATLASZPT7TEV_ptZ.dat
FK: FK_ATLASZPT7TEV_TOT.dat
OP: RATIO" > FK_ATLASZPT7TEV-COMPOUND.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
exit 0
