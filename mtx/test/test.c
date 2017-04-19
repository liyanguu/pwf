#include <stdio.h>
#include "mtx.h"

int main() {
	Spm a, p;
	int i;

	a = makesp(20);
	p = a;
	for (i = 0; i < 20; i++) {
		if (i > 0)
			addacol(p, 1.0, i-1);
		addacol(p, 2.0, i);
		if (i < 19)
			addacol(p, -1.0, i+1);
		p = p->nextrow;
	}
	printsp(a);
}

/* 
#define N 20
Size sn = N;

void asub(Elm *, Elm *);
void atsub(Elm *, Elm *);
void sparse(Elm *b, Size n, void (*asub)(Elm *, Elm *),
	void (*atsub)(Elm*, Elm*), Elm *x, Elm *rsq);
*/
/* test sparse */
int main() {
	Elm b[N], x[N], bcmp[N], rsq;
	Size i;

	for (i = 0; i < N; i++) {
		b[i] = 1.0;
		x[i] = .0;
	}
	b[0] = 3.0;
	b[N-1] = -1.0;
	sparse(b, N, asub, atsub, x, &rsq);
	printf("rsq: %7.5f\n", rsq);
	for (i = 0; i < N; i++)
		printf("%7.5f\n", x[i]);
	asub(x, bcmp);
	for (i = 0; i < N; i++)
		printf("%7.5f  %7.5f\n", bcmp[i], b[i]);
	return 0;
}

void asub(Elm *x, Elm *y) {
	Size i;

	y[0] = x[0] + 2.0 * x[1];
	y[sn-1] = -2.0 * x[sn-2] + x[sn-1];
	for (i = 1; i < sn-1; i++)
		y[i] = -2.0 * x[i-1] + x[i] + 2.0 * x[i+1];
}

void atsub(Elm *x, Elm *y) {
	Size i;

	y[0] = x[0] - 2.0 * x[1];
	y[sn-1] = 2.0 * x[sn-2] + x[sn-1];
	for (i = 1; i < sn-1; i++)
		y[i] = 2.0 * x[i-1] + x[i] - 2.0 * x[i+1];
}
	
/* 
int test_mtx_01(void) {
	Elm a[1000];
	Elm b[ES];
	Elm y[ES], x[ES];
	Mtx m, inv_m;
	int i, j;
	
	for (i = 0; i < 1000 && scanf("%lf", &a[i]) == 1; i++)
		;
	m = mtxalloc(ES, ES);
	inv_m = mtxalloc(ES, ES);
	mtxread(m, a, ES);
	mtxprint("Ybus", m);
	if (lufact(m)) {
		mtxprint("LU", m);
		for (i = 0; i < ES; ++i) {
			for (j = 0; j < ES; ++j)
				b[j] = 0;
			b[i] = 1;
			forward(m, b, y);
			back(m, y, x);
			readcol(inv_m, i, x);
		}	
		mtxprint("Inverse", inv_m);
	}
	return 0;
}
*/
