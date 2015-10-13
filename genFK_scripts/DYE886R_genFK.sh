
if [ ! -e "FK_DYE886R_D.dat" ] 
then
	exit 1;
fi

if [ ! -e "FK_DYE886R_P.dat" ] 
then
	exit 1;
fi

echo "# COMPOUND FK
FK: FK_DYE886R_D.dat
FK: FK_DYE886R_P.dat
OP: RATIO" > FK_DYE886R-COMPOUND.dat

# Cleanup

exit 0
