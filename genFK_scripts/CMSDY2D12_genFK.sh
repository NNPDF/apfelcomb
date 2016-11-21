# Run merge
FKmerge FK_CMSDY2D12_BIN1.dat FK_CMSDY2D12_BIN2.dat FK_CMSDY2D12_BIN3.dat FK_CMSDY2D12_BIN4.dat FK_CMSDY2D12_BIN5.dat FK_CMSDY2D12_BIN6.dat > FK_CMSDY2D12.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup 
rm FK_CMSDY2D12_BIN*
exit 0
