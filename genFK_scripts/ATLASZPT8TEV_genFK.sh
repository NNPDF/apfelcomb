# Increment tables in BIN 4 and BIN6 by 10 dummy points
FKincrement FK_ATLASZPT8TEV-MLLBIN4_ptZ.dat 10 > FK_ATLASZPT8TEV-MLLBIN4_ptZ.tmp.dat
FKincrement FK_ATLASZPT8TEV-MLLBIN6_ptZ.dat 10 > FK_ATLASZPT8TEV-MLLBIN6_ptZ.tmp.dat

## Exclude bin 5 from MLL distribution so that it does not overlap with YDIST!
## Unnormalised distributions
FKmerge FK_ATLASZPT8TEV-MLLBIN1_ptZ.dat FK_ATLASZPT8TEV-MLLBIN2_ptZ.dat FK_ATLASZPT8TEV-MLLBIN3_ptZ.dat FK_ATLASZPT8TEV-MLLBIN4_ptZ.tmp.dat FK_ATLASZPT8TEV-MLLBIN6_ptZ.tmp.dat > FK_ATLASZPT8TEVMDIST.dat

## Merge 6 bins of YLL distribution so that it does not overlap with YDIST!
## Unnormalized distributions
FKmerge FK_ATLASZPT8TEV-YLLBIN1_ptZ.dat FK_ATLASZPT8TEV-YLLBIN2_ptZ.dat FK_ATLASZPT8TEV-YLLBIN3_ptZ.dat FK_ATLASZPT8TEV-YLLBIN4_ptZ.dat FK_ATLASZPT8TEV-YLLBIN5_ptZ.dat FK_ATLASZPT8TEV-YLLBIN6_ptZ.dat > FK_ATLASZPT8TEVYDIST.dat

RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
exit 0
