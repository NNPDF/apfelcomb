theoryID=64

for i in `seq 27 30`; do
    ./apfel_comb dis $i $theoryID
done

for i in `seq 4 23`; do
    ./apfel_comb dyp $i $theoryID
done
