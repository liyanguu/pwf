sources = $(wildcard *.c)
objs = $(patsubst %.c, %.o, $(sources)) simplex/mtx.o
cc = gcc
cflags = -Wall

include depends

pwf: $(objs)
	$(cc) $(cflags) -o $@ $^ -lm 

depends: $(sources)
	gcc -M $^ >$@


simplex/mtx.o: simplex/mtx.c simplex/mtx.h msg.h
