#include "mtx.h"

int main() {
	Mtx a = mtxalloc(3, 3);
	Elm v[] = { 0, 3, 1, 1, 2, -2, 2, 5, 4 };
	Mtx b;

	mtxread(a, v);
	mtxprint("a", a);
	b = mtxdup(a);
	mtxprint("b", b);
	mtxfree(a);
	mtxfree(b);
	return 0;
}
