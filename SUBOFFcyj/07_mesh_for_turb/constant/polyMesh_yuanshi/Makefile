grid.e:grid.o
	ifort -o grid.e grid.o -lpthread -shared-intel
grid.o:grid.f90
	ifort -c -r8 -O3 -mcmodel=large -mtune=pentium4 -assume buffered_io grid.f90