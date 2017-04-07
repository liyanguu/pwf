SRCS = main.c misc.c
OBJS := $(SRCS:%.c=%.o)
DEPS := $(SRCS:%.c=%.d)
CC = gcc
LOADER = gfortran
PWF_PATH = ../pwf
CFLAGS := -Wall -O3
CPPFLAGS := -I$(PWF_PATH)
LDFLAGS := -L$(PWF_PATH)
LA_LIB = lapacke cblas lapack refblas 
KLU_LIB = csparse klu amd colamd btf
LIBS := dl pwf m $(LA_LIB) $(KLU_LIB)
EXES = minq

.PHONY : all pwflib

all : $(EXES)

$(EXES) : $(OBJS) pwflib  
	$(LOADER) -o $@ $(filter-out pwflib, $^) $(LDFLAGS) $(addprefix -l,$(LIBS))

pwflib :
	$(MAKE) -C $(PWF_PATH) lib

%.d : %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< >$@.$$$$; \
	sed 's,\($*\.o\)[ :]*,\1 $@ : ,g' <$@.$$$$ >$@; \
	rm -f $@.$$$$

-include $(DEPS)

.PHONY : clean veryclean
clean :
	rm -rf $(OBJS) $(DEPS)
veryclean : clean
	rm -rf $(EXES)
