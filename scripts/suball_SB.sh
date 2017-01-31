for (( c=0; c<=167; c++ ))
do  
   echo $c
   echo "cd ${PWD}
time ./appl_comb ${c} ${1}" > appl_${c}_${1}.run
   qsub -q generic -l walltime=24:00:00 appl_${c}_${1}.run
done

for (( c=1; c<=43; c++ ))
do
   echo "cd ${PWD}
time ./dis_comb ${c} ${1}" > dis_${c}_${1}.run
   qsub -q generic -l walltime=12:00:00 dis_${c}_${1}.run
done

for (( c=1; c<=23; c++ ))
do
   echo "cd ${PWD}
time ./ftdy_comb ${c} ${1}" > ftdy_${c}_${1}.run
   qsub -q generic -l walltime=12:00:00 ftdy_${c}_${1}.run
done
