# merge FK tables
FKmerge FK_LHCBWMU1FB_WP.dat FK_LHCBWMU1FB_WM.dat > FK_LHCBWMU1FB.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_LHCBWMU1FB_W*

exit 0
