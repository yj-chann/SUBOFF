cylinder.e:cylinder.o
	ifort -o cylinder.e cylinder.o -lpthread -shared-intel
cylinder.o:cylinder.F90
	ifort -c -r8 -O3 -mcmodel=large -mtune=pentium4 -assume buffered_io cylinder.F90