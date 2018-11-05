### Bash script for running extractClusterToPDB.py for each state and running 
### SAXS calculations (using Crysol) for all extracted PDB files for each state on SGE cluster. 
### Required softwares: Crysol
### @Chuankai Zhao, czhao37@illinois.edu

#!/bin/bash

path=$(pwd)

for state in `seq 0 499`
do

cat > PBS_STATE_${state} << EOF
#$ -S /bin/bash    # Set shell to run job
#$ -q all.q        # Choose queue to run job in
#$ -pe orte 1      # Request one processor from the orte parallel env.
#$ -o ${path}/state_${state}.log

cd ${path}
mpirun -np 1 python extractClusterToPDB.py -m MSM25.pkl -c clustering_tica.pkl -ns 100 -np ${state} -s 10

mkdir ${path}/SAXS_STATE_${state}
cd ${path}/SAXS_STATE_${state}
mv ../state${state}_*.pdb .
mpirun -np 1 crysol state${state}_*.pdb -lm 40 -fb 18 -sm 0.5 -ns 51
EOF

done
