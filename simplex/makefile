objs = simp.o ../mtx.o ../msg.o 

test_simp: $(objs) test_simp.o
	gcc $(objs) test_simp.o

test_inverse: mtx.o ../msg.o test_inverse.o 
	gcc -o test_inverse mtx.o ../msg.o test_inverse.o 

test_ludcmp: mtx.o ../msg.o test_ludcmp.o
	gcc mtx.o ../msg.o test_ludcmp.o

test_mtxdup: mtx.o ../msg.o test_mtxdup.o
	gcc -o test_mtxdup mtx.o ../msg.o test_mtxdup.o

test_simp.o: test_simp.c simp.h
test_ludcmp.o: test_ludcmp.c
test_inverse.o: test_inverse.c
test_mtxdup.o: test_mtxdup.c
simp.o: simp.c simp.h
	gcc -c -I.. simp.c
../mtx.o: ../mtx.c
../msg.o: ../msg.c
clear:
	rm $(objs)
