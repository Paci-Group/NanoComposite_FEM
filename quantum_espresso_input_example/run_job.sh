#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0 # memory; default unit is megabytes
#SBATCH --time=0-23:55 # time (DD-HH:MM)
#SBATCH --job-name=Ag8_MgO # Name of job in queue
#SBATCH --account=XXXXXXX
module load nixpkgs/16.09  intel/2018.3  openmpi/3.1.4
module load quantumespresso/6.4


export OMP_NUM_THREADS=1

echo "Starting run at: `date`"

# Electronic Minimization #
mpiexec cp.x < 1_emin_gram.in > 1_emin_gram.out
sleep 10

mpiexec cp.x < 2_emin_damp.in > 2_emin_damp.out
sleep 10


# Ionic and Cell Minimization #
mpiexec cp.x < 3_geom_sd.in > 3_geom_sd.out
sleep 10

mpiexec cp.x < 4_geom_damp.in > 4_geom_damp.out
sleep 10

mpiexec cp.x < 5_geom_cell.in > 5_geom_cell.out
sleep 10

mpiexec cp.x < 6_geom_nstepe.in > 6_geom_nstepe.out
sleep 10

# Dielectric Calculation #

mpiexec cp.x < 7_zero_field_pol.in > 7_zero_field_pol.out
sleep 10

mpiexec cp.x < 8_efield_electron_pol.in > 8_efield_electron_pol.out
sleep 10

mpiexec cp.x < 9_efield_ion_pol.in > 9_efield_ion_pol.out

echo "Program finished with exit code $? at: `date`"

