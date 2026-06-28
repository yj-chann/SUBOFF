#!/bin/bash
#SBATCH -p vip_21
#SBATCH -N 6
#SBATCH -n 384
#SBATCH -J SUB_09_fileGen
#source /es01/paratera/parasoft/module.sh
module purge
module load mpi/intel/19.3.0
mpirun -np 384 ./suboff.zf 
