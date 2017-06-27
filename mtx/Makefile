SRCS = $(wildcard *.c)
OBJS := $(SRCS:.c=.o)
DEPS := $(SRCS:.c=.d)
CC = gcc
CPPFLAGS =
CFLAGS = -Wall -O
RM-F = rm -rf
MTXLIB = libmtx.a

.PHONY : lib clean veryclean

lib : $(MTXLIB)

$(MTXLIB) : $(OBJS)
	ar -rvs $@ $?

%.d : %.c
	$(CC) -MM $(CPPFLAGS) $< >$@

-include $(DEPS)

clean :
	$(RM-F) *.o *.d
veryclean : clean
	$(RM-F) $(MTXLIB)
