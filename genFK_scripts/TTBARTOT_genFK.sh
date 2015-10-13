# Run merge
mergeFK FK_TTBARTOT_7TEV.dat FK_TTBARTOT_7TEV.dat FK_TTBARTOT_7TEV.dat FK_TTBARTOT_7TEV.dat FK_TTBARTOT_8TEV.dat FK_TTBARTOT_8TEV.dat > FK_TTBARTOT_UN.dat 

# Check return value
RETCODE=$?; if [[ $RETCODE != 0 ]]; then exit $RETCODE; fi

# Separate into header and body
awk '1;/{FastKernel/{exit}' FK_TTBARTOT_UN.dat > HEADER.dat
awk 'f;/{FastKernel/{f=1}' FK_TTBARTOT_UN.dat > BODY.dat

# Bin width and unit correction (*14000/1000)
awk ' CONVFMT="%.16g"; OFMT="%.16g"; { for(i = 4; i <= NF; i++) { $i=$i*14.0; }; print; } ' BODY.dat > BODY_NM.dat
cat HEADER.dat BODY_NM.dat > FK_TTBARTOT.dat

rm HEADER.dat BODY.dat BODY_NM.dat FK_TTBARTOT_UN.dat FK_TTBARTOT_7TEV.dat FK_TTBARTOT_8TEV.dat

exit 0
