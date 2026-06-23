mkdir.e:mkdir.o
	ifort -o mkdir.e mkdir.o -lpthread -shared-intel
mkdir.o:mkdir.F90
	ifort -c -r8 -O3 -mcmodel=large -mtune=pentium4 -assume buffered_io mkdir.F90