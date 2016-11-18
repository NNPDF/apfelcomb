# Run merge
FKmerge FK_CMSDY2D11_BIN1.dat FK_CMSDY2D11_BIN2.dat FK_CMSDY2D11_BIN3.dat FK_CMSDY2D11_BIN4.dat FK_CMSDY2D11_BIN5.dat FK_CMSDY2D11_BIN6.dat > FK_CMSDY2D11.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_CMSDY2D11_BIN*
exit 0
