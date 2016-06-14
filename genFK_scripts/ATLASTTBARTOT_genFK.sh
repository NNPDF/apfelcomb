# merge FK tables
FKmerge FK_ATLASTTBARTOT7TEV.dat ATLASTTBARTOT8TEV.dat ATLASTTBARTOT13TEV.dat > FK_ATLASTTBARTOT.dat

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Cleanup
rm FK_ATLASTTBARTOT7TEV.dat ATLASTTBARTOT8TEV.dat ATLASTTBARTOT13TEV.dat

exit 0
