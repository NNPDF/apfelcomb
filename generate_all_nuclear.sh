theory_ID=65

##DIS_NC
for i in {96..166}
do
   ./apfel_comb dis $i $theory_ID
done
 
##DIS_CC
./apfel_comb dis 5 $theory_ID
./apfel_comb dis 6 $theory_ID
./apfel_comb dis 23 $theory_ID
./apfel_comb dis 24 $theory_ID
 
##APP
 for i in {2062..2080}
 do
    ./apfel_comb app $i $theory_ID
 done
 
for i in {2081..2090}
do
   ./apfel_comb app $i $theory_ID
done

##Positivity
##F2
for i in {54..58}
do
   ./apfel_comb dis $i $theory_ID
done
#DY
for i in {4..7}
do
   ./apfel_comb dyp $i $theory_ID
done
for i in {16..23}
do
   ./apfel_comb dyp $i $theory_ID
done

./apfel_comb dis 1 $theory_ID #SLACD (Remember to set DIS_F2P)
./apfel_comb dis 19 $theory_ID #NMCPD_P
./apfel_comb dis 25 $theory_ID #SLACD (Remember to set DIS_F2P)

./apfel_comb dyp 0 $theory_ID #DYE605
./apfel_comb dyp 1 $theory_ID #DYE886R_P