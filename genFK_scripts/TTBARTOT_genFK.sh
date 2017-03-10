# Run merge
FKmerge FK_TTBARTOT_7TEV.dat FK_TTBARTOT_7TEV.dat FK_TTBARTOT_7TEV.dat FK_TTBARTOT_7TEV.dat FK_TTBARTOT_8TEV.dat FK_TTBARTOT_8TEV.dat > FK_TTBARTOT.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

exit 0
