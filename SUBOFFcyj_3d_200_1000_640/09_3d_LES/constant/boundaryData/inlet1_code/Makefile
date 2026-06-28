suboff_main.e:suboff_main.o
	ifort -o suboff_main.e suboff_main.o -lpthread -shared-intel
suboff_main.o:suboff_main.F90
	ifort -c -r8 -O3 -mcmodel=large -mtune=pentium4 -assume buffered_io suboff_main.F90