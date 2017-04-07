CC := gcc
LOADER := gfortran
CFLAGS := -Wall -O -I..
LIBS := m lapacke cblas lapack refblas
EXES := test_selectrow test_inverse test_mtxinv
vpath %.h ..

.PHONY : all
all : $(EXES)

$(EXES) : % : %.o ../mtx.o ../msg.o
	$(LOADER) -o $@ $^ $(addprefix -l, $(LIBS))

%.o : %.c mtx.h msg.h
	$(CC) -c $(CFLAGS) -o $@ $<

.PHONY : clean
clean :
	rm -rf *.o
