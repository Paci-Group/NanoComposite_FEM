#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00-02:55 # time (DD-HH:MM)
#SBATCH --job-name=ag12 # Name of job in queue
#SBATCH --account=XXXXX
module load CCEnv
module load StdEnv/2020  intel/2020.1.217  openmpi/4.0.3
module load siesta/4.1-b4

mpirun siesta < file.fdf > file.out
