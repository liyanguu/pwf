#include <stdio.h> 
#include <stdlib.h> 
typedef double Elm;
typedef int Size;
typedef struct Mtx {
	Elm *v;
	Size nrow, ncol;
} Mtx;
typedef struct Dia {
	Elm *v;
	Size n;
} Dia;

int mtxalloc(Mtx *m, Size nr, Size nc);
Size pos(Mtx m, Size i, Size j);
Elm mtxget(Mtx m, Size i, Size j);
int mtxset(Elm e, Mtx m, Size i, Size j);
void mtxprn(char *s, Mtx m);
void lufact(Mtx a);
void mtxread(Mtx m, Elm *v);
void fwsub(Mtx a, Elm *v, Elm *y);

main() {
	Mtx a, yv;
	Elm data[9] = { 1, 1, 1, 2, 4, 1, -3, 1, -2 };
	Elm b[3] = { 4, 3, -13 };
	Elm y[3], x[3];
	
	mtxalloc(&a, 3, 3);
	mtxread(a, data);
	mtxprn("A", a);
	lufact(a);
	mtxprn("LU", a);
	fwsub(a, b, y);
	mtxalloc(&yv, 3, 1);
	mtxread(yv, y);
	mtxprn("y", yv);
}
	
int mtxalloc(Mtx *m, Size nr, Size nc) {
	if ((m->v = (Elm *)malloc(sizeof(Elm) * nr * nc)) == NULL)
		return -1;
	m->nrow = nr;
	m->ncol = nc;
	return 0;
}

Size pos(Mtx m, Size i, Size j) {
	if (i < 0 || i >= m.nrow || j < 0 || j >= m.ncol)
		return -1;
	return i * m.ncol + j;
}

Elm mtxget(Mtx m, Size i, Size j) {
	Size k;

	if ((k = pos(m, i, j)) < 0)
		return 0.0;
	return m.v[k];
} 

int mtxset(Elm e, Mtx m, Size i, Size j) {
	Size k;

	if ((k = pos(m, i, j)) < 0)
		return -1;
	m.v[k] = e;
	return 0;
}

void mtxread(Mtx m, Elm *v) {
	Size i, j;
	
	for (i = 0; i < m.nrow; ++i)
		for (j = 0; j < m.ncol; ++j)
			mtxset(*v++, m, i, j);
}

void mtxprn(char *s, Mtx m) {
	Size i, j;

	printf("%s\n", s);
	for (i = 0; i < m.nrow; ++i) {
		for (j = 0; j < m.ncol; ++j)
			printf("%8.3f ", (double) mtxget(m, i, j));
		putchar('\n');
	}
	putchar('\n');
}

void mtxmul(Mtx m, Mtx m1, Mtx m2) {
	
}

void lufact(Mtx a) {
	Elm m, x, y;
	Size i, k, j;

	for (i = 0; i < a.nrow - 1; ++i)
		for (k = i+1; k < a.nrow; ++k) {
			m = mtxget(a, k, i) / mtxget(a, i, i);
			mtxset(m, a, k, i);
			for (j = i+1; j < a.ncol; ++j) {
				x = mtxget(a, k, j);
				y = mtxget(a, i, j); 
				mtxset(x - m * y, a, k , j); 
			}
		}
}

void fwsub(Mtx a, Elm *v, Elm *y) {
	Size i, j;

	for (i = 0; i < a.nrow; ++i)
		for (y[i] = v[i], j = 0; j < i; ++j)
			y[i] -= mtxget(a, i, j);
}
