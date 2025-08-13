#!/bin/bash

# Slurm job options (name, compute nodes, job time)

#SBATCH --job-name=Blaze
#SBATCH --time=0:20:0

#SBATCH --cpus-per-task=1


#SBATCH --account=mdisspt-s2259894
#SBATCH --partition=standard
#SBATCH --qos=standard 

module load intel-20.4/compilers
module load intel-20.4/mpi
module load intel-20.4/cmkl

 # Change to the submission directory
cd $SLURM_SUBMIT_DIR

srun -n $1 ../bin/Blaze --nsteps 100 --nbays 140 --nfloors 5 --udl -2.2e5  --elem_type LinearElastic
