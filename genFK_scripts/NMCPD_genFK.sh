if [ ! -e "FK_NMCPD_D.dat" ] 
then
	exit 1;
fi

if [ ! -e "FK_NMCPD_P.dat" ] 
then
	exit 1;
fi

echo "# COMPOUND FK
FK: FK_NMCPD_D.dat
FK: FK_NMCPD_P.dat
OP: RATIO" > FK_NMCPD-COMPOUND.dat

# Cleanup

exit 0
