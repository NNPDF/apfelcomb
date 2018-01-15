# This script is the submission script for the stoomboot cluster
echo "cd /project/theorie/nhartlan/apfelcomb
./apfel_comb $1 $2 $3" > apc_${1}_${2}_${3}.run
qsub -q generic -l walltime=24:00:00 apc_${1}_${2}_${3}.run
