# 用户可配置选项
FC          = mpiifort
FLINKER     = mpiifort
OPTFLAGS    = -w -mcmodel=large -r8 -O3
WARNFLAGS   = -w

FFLAGS = $(OPTFLAGS) $(WARNFLAGS)
EXECS = suboff.zf

# 明确列出所有源文件（按依赖顺序可选）
# SRCS = suboff_main.F90 \
#        grid.F90 \
#        ini_up_points.F90 \
#        code_main.F90 \
#        projection.F90

SRCS = suboff_main.F90 \
   	   grid.F90 \
	   ini_up_points.F90 \
	   code_main.F90
OBJS = $(SRCS:.F90=.o)

# 默认目标
default: $(EXECS)

# 链接生成可执行文件
$(EXECS): $(OBJS)
	$(FLINKER) $(FFLAGS) -o $@ $^

# 编译规则
%.o: %.F90
	$(FC) $(FFLAGS) -c $<

# 清理
clean:
	/bin/rm -f *.o $(EXECS) *.out

# 头文件依赖（当 head.fi 或 parameter.h 修改时，强制重新编译所有源文件）
head.fi parameter.h common.fi param.fi:
	@touch $(SRCS)