#!/bin/bash
#SBATCH -p vip_21
#SBATCH -N 1
#SBATCH -n 50
#SBATCH -J SUB_06
module purge
module load mpi/intel/19.3.0
cd /public1/home/scb7552/cyj/SUBOFF/SUBOFFMeshConstruct/06_tg/
mpirun -np 50 ./tg_zf.out




