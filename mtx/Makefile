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
	@set -e; $(RM-F) $@; \
	$(CC) -MM $(CPPFLAGS) $< >$@.$$$$; \
	sed 's/\($*\.o\)[ :]*/\1 $@ : /g' $@.$$$$ >$@; \
	$(RM-F) $@.$$$$; \
	echo "generated $@"

-include $(DEPS)

clean :
	$(RM-F) *.o *.d
veryclean : clean
	$(RM-F) $(MTXLIB)
