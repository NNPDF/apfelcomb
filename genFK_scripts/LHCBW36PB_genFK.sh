# merge FK tables
mergeFK FK_LHCBW36PB_WP.dat FK_LHCBW36PB_WM.dat > FK_LHCBW36PB.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_LHCBW36PB_W*

exit 0
