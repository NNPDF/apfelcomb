source ~/.zshrc

theoryID=64

#EIC
for id in `seq 84 119`; do 
    ./apfel_comb dis $id $theoryID
done