file_out_manualDecomposition.e:file_out_manualDecomposition.o
	ifort -o file_out_manualDecomposition.e file_out_manualDecomposition.o -lpthread -shared-intel
file_out_manualDecomposition.o:file_out_manualDecomposition.f90
	ifort -c -r8 -O3 -mcmodel=large -mtune=pentium4 -assume buffered_io file_out_manualDecomposition.f90