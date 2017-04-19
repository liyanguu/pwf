SRCS = main.c misc.c
OBJS := $(SRCS:%.c=%.o)
DEPS := $(SRCS:%.c=%.d)
CC = gcc
LOADER = gfortran
PWF_PATH = ../pwf
MTX_PATH = ../mtx
CFLAGS := -Wall -O3
CPPFLAGS := -I$(PWF_PATH) -I$(MTX_PATH)
LDFLAGS := -L$(PWF_PATH) -L$(MTX_PATH)
LA_LIB = lapacke cblas lapack refblas 
KLU_LIB = csparse klu amd colamd btf
LIBS := dl pwf mtx m $(LA_LIB) $(KLU_LIB)
EXES = minq

.PHONY : all

all : $(EXES)

$(EXES) : $(OBJS)
	$(LOADER) -o $@ $(OBJS) $(LDFLAGS) $(addprefix -l,$(LIBS))

%.d : %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< | \
	sed 's,\($*\.o\)[ :]*,\1 $@ : ,g' >$@

-include $(DEPS)

.PHONY : clean veryclean
clean :
	rm -rf $(OBJS)
veryclean : clean
	rm -rf $(EXES) $(DEPS)
