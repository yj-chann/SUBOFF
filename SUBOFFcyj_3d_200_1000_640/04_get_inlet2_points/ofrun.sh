#!/bin/bash
#SBATCH -p vip_21
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -J SUB_04_inlet2_runOF


module purge
#module load mpi/intel/19.3.0
source /public1/home/scb7552/soft/OpenFOAM-v2312/env.sh
# blockMesh

writeMeshObj -constant 
rm mesh*.obj patch_OUTLET_constant.obj patch_HULL_constant.obj