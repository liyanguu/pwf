CC := gcc
LOADER := gfortran
CFLAGS := -Wall -O
CPPFLAGS := -I.. 
LIBS := m lapacke cblas lapack refblas
## EXES := test_sp test_selectrow test_inverse test_mtxinv
EXES := test_colspm 
vpath %.h .. 

.PHONY : all
all : $(EXES)

$(EXES) : % : %.o ../colspm.o ../mtx.o ../msg.o
	$(LOADER) -o $@ $^ $(addprefix -l, $(LIBS))

%.o : %.c mtx.h msg.h colspm.h
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

.PHONY : clean
clean :
	rm -rf *.o
