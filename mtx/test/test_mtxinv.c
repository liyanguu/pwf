#include <stdio.h>
#include "mtx.h"

int main() {
	Elm va[] = { 8, 3, 1, 1, 2, -2, 2, 5, 4 };
	Mtx a = mtxalloc(3, 3);
	Mtx b = mtxalloc(3, 3);
	Mtx inva;

	mtxread(a, va);
	mtxprint("a", a);
	inva = mtxinv(a);
	mtxprint("a ** -1", inva);
	mtxmul(b, a, inva);
	mtxprint("ans", b);
}
