CC := gcc
LOADER := gfortran
CFLAGS := -Wall -O
CPPFLAGS := -I.. -I../..
LIBS := m lapacke cblas lapack refblas
EXES := test_sp test_selectrow test_inverse test_mtxinv
vpath %.h .. ../..

.PHONY : all
all : $(EXES)

$(EXES) : % : %.o ../spm.o ../mtx.o ../../msg.o
	$(LOADER) -o $@ $^ $(addprefix -l, $(LIBS))

%.o : %.c spm.h mtx.h msg.h
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

.PHONY : clean
clean :
	rm -rf *.o
