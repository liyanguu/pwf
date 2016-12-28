#include "mtx.h"

/* test function inverse() */
int main() {
	Mtx a = mtxalloc(3, 3);
	Mtx inv;
	Mtx ans;
	Elm v[] = { 0, 3, 1, 1, 2, -2, 2, 5, 4 };

	mtxread(a, v);
	inv = inverse(a);
	mtxprint("inverse", inv);
	ans = mtxmul(a, inv);
	mtxprint("ans", ans);
}
