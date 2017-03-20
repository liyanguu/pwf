SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRCS))
CC = gcc
CFLAGS = -Wall -O3
LIBS := m cblas lapacke blas lapack

pwf : $(OBJS)
	$(CC) -o $@ $^ $(addprefix -l, $(LIBS))

depends : $(SRCS)
	$(CC) -MM $^ >$@

include depends

.PHONY : clean
clean : 
	rm -rf $(OBJS) depends pwf
