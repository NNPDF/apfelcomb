#HADRONIC
for id in `seq 1813 1817`; do
./apfel_comb app $id 64
done

#DIS
for id in `seq 44 67`; do
./apfel_comb dis $id 64
done

#DY POS
for id in `seq 4 6`; do
./apfel_comb dyp $id 64
done

#DIS POS
for id in `seq 27 30`; do
./apfel_comb dis $id 64
done