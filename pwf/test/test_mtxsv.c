#include <stdio.h>
#include "mtx.h"

int main() {
	Elm va[] = { 0, 3, 1, 1, 2, -2, 2, 5, 4 };
	Elm vb[] = { 1, 7, -1 };
	Mtx a = mtxalloc(3, 3);
	Mtx acopy;
	Mtx b = mtxalloc(3, 1);
	Mtx c = mtxalloc(3, 3);
	int flag;

	mtxread(a, va);
	acopy = mtxdup(a);
	mtxprint("a", a);
	mtxread(b, vb);
	flag = mtxsv(a, b);
	mtxprint("LU", a);
	mtxprint("solve", b);
	printf("flag = %d\n", flag);
	mtxmul(c, acopy, b);
	mtxprint("c", c);
}
