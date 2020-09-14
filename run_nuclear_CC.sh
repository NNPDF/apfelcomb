source ~/.zshrc

theoryID=64

#DIS
for id in `seq 44 47`; do
    ./apfel_comb dis $id $theoryID
done