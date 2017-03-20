SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(sources))
CC = GCC
CFLAGS = -Wall

include depends

pwf : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm 

$(OBJS) : $(SRCS)

depends : $(sources)
	$(CC) -M $^ >$@

.PHONY clean
clean : 
	rm -rf $(OBJS) depends


