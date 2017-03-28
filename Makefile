SRCS = main.c pwf.c cdf.c pf.c comp.c msg.c 
OBJS = $(patsubst %.c, %.o, $(SRCS))
CC = gcc
CFLAGS = -Wall -O3
KLU_PATH = klusolve-code-45
LIBS = m klu csparse amd btf colamd

vpath %.h $(KLU_PATH)/Include

pwf : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -L$(KLU_PATH)/Lib $(addprefix -l, $(LIBS))

main.o pwf.o cdf.o pf.o : pwf.h msg.h
pwf.o cdf.o : comp.h
pf.o : pf.c comp.h cs.h klu.h
	$(CC) -c $(CFLAGS) -o $@ $< -I$(KLU_PATH)/Include

depends : $(SRCS)
	$(CC) -MM $^ >$@ -I$(KLU_PATH)/Include

-include depends

.PHONY : clean
clean : 
	rm -rf $(OBJS)
