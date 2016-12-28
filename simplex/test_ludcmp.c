#include <stdio.h>
#include "mtx.h"

int main() {
	Elm v[] = { 0, 3, 1, 1, 2, -2, 2, 5, 4 };
	Elm b[] = { 1, 7, -1 };
	Mtx a = mtxalloc(3, 3);
	Elm det;
	Size indx[100];
	int i;

	mtxread(a, v);
	ludcmp(a, indx, &det);
	mtxprint("LU", a);
	lubksb(a, indx, b);
	for (i = 0; i < 3; i++)
		printf("%f ", b[i]);
}
