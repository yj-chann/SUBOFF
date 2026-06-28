#!/bin/bash
rm -rf *.o
rm -rf *.zf
make

for i in {0..383}; do
    mkdir -p  "../processor$i/constant/polyMesh"  "../processor$i/0" 
done

sbatch ./ofrun.sh
