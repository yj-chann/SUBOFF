
FC = mpiifort
OPT = -r8 -O3

EXEC = tg_zf.out

FILELIST = tg_parallel.F90 

#FFTLIB = -DMKL_DOUBLE \
#-I"/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl"/include \
#-I"/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl"/include/fftw
#
#MKLLIB = -I"/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl"/include \
#-I"/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl"/include/fftw\
#-L/es01/paratera/sce4049/zf/fftw2xf/lib/intel64 -lfftw2xf_double_intel \
#-L"/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin" \
#"/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin"/libmkl_intel_lp64.so \
#-Wl,--start-group "/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin"/libmkl_intel_thread.so \
#"/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin"/libmkl_core.so \
#-Wl,--end-group -L"/es01/software/intel2019/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin" \
#-liomp5 -lpthread -lm 

#in fact, only two paths are altered.
FFTLIB = -DMKL_DOUBLE \
-I"/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl"/include \
-I"/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl"/include/fftw

MKLLIB = -I"/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl"/include \
-I"/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl"/include/fftw\
-L/public1/home/scb7552/njy21/SUBOFF/fftw2xf/lib/intel64 -lfftw2xf_double_intel \
-L"/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl/lib/intel64_lin" \
"/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl/lib/intel64_lin"/libmkl_intel_lp64.so \
-Wl,--start-group "/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl/lib/intel64_lin"/libmkl_intel_thread.so \
"/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl/lib/intel64_lin"/libmkl_core.so \
-Wl,--end-group -L"/public1/soft/intel/2017/compilers_and_libraries_2017.7.259/linux/mkl/lib/intel64_lin" \
-liomp5 -lpthread -lm 



OBJS = $(FILELIST:.F90=.o)

%.o:%.F90
	$(FC) $(OPT) $(FFTLIB) -c  $<


default: $(OBJS)
	$(FC) -r8 -O3 $(FFTLIB) $(OBJS) $(MKLLIB) -o $(EXEC)

clean:
	rm -f *.out *.o work.* *.pc *.mod

