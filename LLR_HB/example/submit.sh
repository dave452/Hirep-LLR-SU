#!/bin/bash

#SBATCH --job-name=3.624.256.pre
#SBATCH --account=scw1019
#SBATCH --ntasks=256

 
###
module load mpi
module load compiler/intel/2018/4 mpi/intel/2018/4
bash setup_replicas.sh -r %n -A pre.dat
mpirun ./llr_hb
