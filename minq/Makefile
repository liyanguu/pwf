SRCS = main.c misc.c
OBJS := $(SRCS:.c:=.o)
DEPS := $(SRCS:.c:=.d)
CC = gcc
LOADER = gfortran
CFLAGS = -Wall -O3 -I..
CPPFLAGS = -I..
LA_LIB = lapacke cblas lapack refblas 
KLU_LIB = csparse klu amd colamd btf
LIBS := dl pwf m $(LA_LIB) $(KLU_LIB)

minq: $(OBJS) 
	$(MAKE) -C .. lib
	$(LOADER) -o $@ $^ -L.. $(addprefix -l,$(LIBS))

include $(DEPS)

%.d : %.c
	@set -e ;\
	rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< >$@.$$$$; \
	sed 's,\($*\.o\)[ :]*,\1 $@ : ,g' <$@.$$$$ >$@; \
	rm -f $@.$$$$

.PHONY : clean
clean :
	rm -rf $(OBJS) $(DEPS)
