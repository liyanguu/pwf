SRCS := $(wildcard *.c)
OBJS := $(patsubst %.c, %.o, $(SRCS))
DEPS := $(patsubst %.c, %.d, $(SRCS))
CC = gcc
LOADER = gfortran
CFLAGS = -Wall -O3 
MTX_PATH = ../mtx
CPPFLAGS = -I$(MTX_PATH)
LDFLAGS = -L$(MTX_PATH)
LA_LIB = lapacke cblas lapack refblas tmglib
KLU_LIB = klu csparse amd btf colamd 
LIBS := mtx $(LA_LIB) $(KLU_LIB) m

.PHONY : all lib mtxlib
all : pwf lib

pwf : $(OBJS) mtxlib
	$(LOADER) -o $@ $(OBJS) $(LDFLAGS) $(addprefix -l, $(LIBS))

mtxlib :
	$(MAKE) -C $(MTX_PATH) lib

lib : libpwf.a

libpwf.a : $(filter-out main.o, $(OBJS))
	ar crs $@ $?

depends : $(SRCS)
	$(CC) -MM $(CPPFLAGS) $^ >$@ 

-include depends

.PHONY : veryclean clean
clean : 
	rm -rf $(OBJS)
veryclean : clean
	rm -rf pwf depends libpwf.a
