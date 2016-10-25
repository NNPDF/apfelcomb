for ((i=1;i<=3;i++)); do
	## Merge together the Z total xsec grids N_applgrid times
	FKselfmerge FK_ATLASZPT7TEV-BIN${i}_tot.dat 14 > FK_ATLASZPT7TEV-BIN${i}_tot_MRG.dat

	# Add 12 dummy points below the Z pT cut
	FKincrement FK_ATLASZPT7TEV-BIN${i}_tot_MRG.dat 12 > FK_ATLASZPT7TEV-BIN${i}_tot_MRG_INC.dat
	FKincrement FK_ATLASZPT7TEV-BIN${i}_ptZ.dat 12 > FK_ATLASZPT7TEV-BIN${i}_ptZ_INC.dat
done 

## Merge the total grids
FKmerge FK_ATLASZPT7TEV-BIN1_tot_MRG_INC.dat FK_ATLASZPT7TEV-BIN2_tot_MRG_INC.dat FK_ATLASZPT7TEV-BIN3_tot_MRG_INC.dat > FK_ATLASZPT7TEV_TOT.dat
FKmerge FK_ATLASZPT7TEV-BIN1_ptZ_INC.dat FK_ATLASZPT7TEV-BIN2_ptZ_INC.dat FK_ATLASZPT7TEV-BIN3_ptZ_INC.dat > FK_ATLASZPT7TEV_PTZ.dat

# Remove the temporary grids
rm FK_ATLASZPT7TEV-BIN*

## write the compound file for each bin
echo "# COMPOUND FK
FK: FK_ATLASZPT7TEV_PTZ.dat
FK: FK_ATLASZPT7TEV_TOT.dat
OP: RATIO" > FK_ATLASZPT7TEV-COMPOUND.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
exit 0
