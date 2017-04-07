#include "mtx.h"

/* test function inverse() */
int main() {
	Mtx a = mtxalloc(3, 3);
	Mtx inva;
	Mtx ans = mtxalloc(3, 3);
	Elm v[] = { 9, 3, 1, 1, 2, -2, 2, 5, 4 };
	Elm vs[] = { 5, 3, 1, 3, 2, -2, 1, -2, 4};

	mtxread(a, vs);
	mtxprint("a", a);
	inva = inv(a);
	mtxprint("inverse", inva);
	mtxmul(ans, a, inva);
	mtxprint("ans", ans);

	mtxread(a, vs);
	mtxprint("a", a);
	inva = sinv(a);
	mtxprint("inverse", inva);
	mtxmul(ans, a, inva);
	mtxprint("ans", ans);

	inva = reshape(inva, 5, 6);
	mtxprint("5 x 6", inva);
	inva = reshape(inva, 3, 3);
	mtxprint("3 x 3", inva);
	inva = reshape(inva, 2, 4);
	mtxprint("2 x 4", inva);
	return 0;
}
