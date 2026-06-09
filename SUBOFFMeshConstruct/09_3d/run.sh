# mpirun --use-hwthread-cpus -np 16 redistributePar -decompose -cellDist -parallel -dry-run
cp -r 0.orig 0 && cp -r constant/polyMesh.orig constant/polyMesh
mpirun --use-hwthread-cpus -np 16 redistributePar -decompose -cellDist -parallel -overwrite
mpirun --use-hwthread-cpus -np 16 topoSet -parallel
mpirun --use-hwthread-cpus -np 16 refineMesh -parallel  -overwrite 
for i in {0..15} ; do rm -rf "processor$i/0//polyMesh"  ; done
mpirun --use-hwthread-cpus -np 16 redistributePar -reconstruct -parallel -overwrite -constant # only reconstruct static mesh to case/constant/polyMesh
reconstructPar -time "0.001"
mpirun --use-hwthread-cpus -np 16 redistributePar -reconstruct -parallel -overwrite -constant
foamCleanPolyMesh
for i in {0..15} ; do mv "processor$i" hidden.procs ; done
for i in {0..15} ; do mv  "hidden.procs/processor$i"  . ; done

mpirun --use-hwthread-cpus -np 16 pimpleFoam  -parallel > log_run

for i in {0..15} ; do for t in 0.1 0.2 0.3 0.4 0.5; do rm -rf  "processor$i/${t}"  . ; done ; done
