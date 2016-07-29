# Increment tables in BIN 4 and BIN6 by 10 dummy points
FKincrement FK_ATLASZPT8TEV-MLLBIN4_ptZ.dat 10
FKincrement FK_ATLASZPT8TEV-MLLBIN6_ptZ.dat 10

## ONLy mll distributions so far - exclude bin 5 so not to overlap with YDIST!
FKmerge FK_ATLASZPT8TEV-MLLBIN1_ptZ.dat FK_ATLASZPT8TEV-MLLBIN2_ptZ.dat FK_ATLASZPT8TEV-MLLBIN3_ptZ.dat FK_ATLASZPT8TEV-MLLBIN4_ptZ.dat FK_ATLASZPT8TEV-MLLBIN6_ptZ.dat > FK_ATLASZPT8TEVMDIST.dat
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
exit 0
