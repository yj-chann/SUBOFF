#!/bin/bash
#SBATCH -p vip_21
#SBATCH -N 2
#SBATCH -n 128
#SBATCH -J SUB_05_runOF


module purge
#module load mpi/intel/19.3.0
# source /public1/home/scb7552/soft/OpenFOAM-v2312/env.sh
source /public1/home/scb7552/openform-6-3/OpenFOAM-v2406/etc/rebashrc

# foamCleanPolyMesh
# blockMesh
mpirun -np 128 redistributePar -decompose -parallel -overwrite
mpirun -np 128 simpleFoam  -parallel
# mpirun -np 64 simpleFoam -postProcess -func "turbulenceFields(fields=(R))" -parallel -latestTime
mpirun -np 128 redistributePar -reconstruct -parallel -overwrite -latestTime
# mv "50000/turbulenceProperties:R" "50000/R"

# for i in {0..127} ; do mv "processor$i" hidden.procs ; done
# for i in {0..127} ; do mv "hidden.procs/processor$i" "processor$i" ; done
# postProcess -func "components(U)" -latestTime
# postProcess -func writeCellCentres -latestTime
# simpleFoam -postProcess -func wallShearStress -latestTime
# find processor* -maxdepth 1 -mindepth 1 -type d ! -name "0" ! -name "6000" -regextype posix-extended -regex ".*/[0-9]+(\.[0-9]+)?$" -exec rm -rf {} + 
# 删除processor中间时间目录
