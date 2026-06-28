cylinder.e:cylinder.o
	ifort -o cylinder.e cylinder.o -lpthread -shared-intel
# -shared-intel 生成的可执行文件依赖共享库
cylinder.o:cylinder.F90
	ifort -c -r8 -O3 -mcmodel=large -mtune=pentium4 -assume buffered_io cylinder.F90 
# -r8 将所有的 REAL 变量默认声明为8字节 -O3 开启最高级别的优化（包括循环展开、向量化等）-mcmodel=large 指定代码和数据的内存模型为“large”