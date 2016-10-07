## Merge 6 bins of YLL distribution so that it does not overlap with YDIST!
## Unnormalized distributions
FKmerge FK_ATLASZPT8TEV-YLLBIN1_ptZ.dat FK_ATLASZPT8TEV-YLLBIN2_ptZ.dat FK_ATLASZPT8TEV-YLLBIN3_ptZ.dat FK_ATLASZPT8TEV-YLLBIN4_ptZ.dat FK_ATLASZPT8TEV-YLLBIN5_ptZ.dat FK_ATLASZPT8TEV-YLLBIN6_ptZ.dat > FK_ATLASZPT8TEVYDIST.dat

RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
rm FK_ATLASZPT8TEV-YLLBIN*

exit 0
