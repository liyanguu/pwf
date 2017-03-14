#include "mtx.h"

/* test function inverse() */
int main() {
	Mtx a = mtxalloc(3, 3);
	Mtx inv;
	Mtx ans = mtxalloc(3, 3);
	Elm v[] = { 0, 3, 1, 1, 2, -2, 2, 5, 4 };

	mtxread(a, v);
	mtxprint("a", a);
	inv = inverse(a);
	mtxprint("inverse", inv);
	mtxmul(ans, a, inv);
	mtxprint("ans", ans);

	inv = reshape(inv, 5, 6);
	mtxprint("5 x 6", inv);
	inv = reshape(inv, 3, 3);
	mtxprint("3 x 3", inv);
	inv = reshape(inv, 2, 4);
	mtxprint("2 x 4", inv);
	return 0;
}
