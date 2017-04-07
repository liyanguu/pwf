SRCS = main.c pwf.c cdf.c pf.c comp.c msg.c mtx.c
OBJS = $(patsubst %.c, %.o, $(SRCS))
CC = gcc
LOADER = gfortran
CFLAGS = -Wall -O3 
LA_LIB = lapacke cblas lapack refblas tmglib
KLU_LIB = klu csparse amd btf colamd 
LIBS := $(LA_LIB) $(KLU_LIB) m 

.PHONY : all lib
all : pwf lib

pwf : $(OBJS)
	$(LOADER) -o $@ $^ $(addprefix -l, $(LIBS))

lib : libpwf.a

libpwf.a : $(filter-out main.o, $(OBJS))
	ar crs $@ $?

depends : $(SRCS)
	$(CC) -MM $^ >$@ 

-include depends

.PHONY : veryclean clean
clean : 
	rm -rf $(OBJS)
veryclean : clean
	rm -rf pwf depends libpwf.a
