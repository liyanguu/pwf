#include <stdio.h>
#include "mtx.h"

int main() {
	Elm va[] = { 0, 3, 1, 1, 2, -2, 2, 5, 4 };
	Elm vb[] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
	Mtx a = mtxalloc(3, 3);
	Mtx b = mtxalloc(3, 3);
	Mtx c = mtxalloc(3, 3);
	Mtx d;

	mtxread(a, va);
	mtxread(b, vb);
	mtxprint("a", a);
	mtxprint("b", b);
	d = trans(a);
	mtxmul(c, a, b);
	mtxprint("c", c);
	xmtxmul(c, d, b);
	mtxprint("c", c);
	return 0;
}
