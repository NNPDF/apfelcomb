source ~/.zshrc

theoryID=64

#HADRONIC
for id in `seq 1813 1817`; do
    ./apfel_comb app $id $theoryID
done

#DIS
for id in `seq 44 67`; do
    ./apfel_comb dis $id $theoryID
done

#DY POS
for id in `seq 4 23`; do
    ./apfel_comb dyp $id $theoryID
done

#DIS POS
for id in `seq 27 30`; do
    ./apfel_comb dis $id $theoryID
done