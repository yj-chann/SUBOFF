#!/bin/bash
rm -rf tg_zf.out
rm -rf tg_parrallel.o
make
if [ ! -f tg_zf.out ]; then
    echo "Error: Compilation failed, tg_zf.out not found."
    exit 1
fi
sbatch ofrun.sh 

