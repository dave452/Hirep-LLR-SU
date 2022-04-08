#!/bin/bash

#SBATCH --job-name=su3.paratest.4x20_0	
#SBATCH --account=scw1019
#SBATCH --ntasks=256
#SBATCH --time=0-10:00

 
###
module load mpi
module load compiler/intel/2018/4 mpi/intel/2018/4

mpirun ./llr_hb
