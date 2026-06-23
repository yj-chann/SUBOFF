#!/bin/bash
#SBATCH -p vip_21
#SBATCH -N 2
#SBATCH -n 100
#SBATCH -J SUB_07_tgGen
#source /es01/paratera/parasoft/module.sh
module purge
module load mpi/intel/19.3.0
mpirun -np 100 ./a.zf  
rm -rf a.zf
rm -rf suboff_main.o
