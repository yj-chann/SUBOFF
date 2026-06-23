#!/bin/bash
rm -rf a.zf
rm -rf suboff_main.o
make
if [ ! -f a.zf ]; then
    echo "Error: Compilation failed, tg_zf.out not found."
    exit 1
fi
sbatch ofrun.sh 

