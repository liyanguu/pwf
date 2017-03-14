#include <stdio.h>
#include "mtx.h"

Spm makesp(Size, Size);
void add(Spm, Size, Size, Elm);
void printsp(char *, Spm);
void deletsp(Spm);
void upper(Spm);
void gauss_jordan(Spm, Elm *);
void printv(Elm *, Size);

int main() {
	Spm a;
	Mtx ma;
	int i, j;
	Elm b[] = { 5, 6, 6, 6, 6, 6, 3 };

	a = makesp(7, 7);
	ma = mtxalloc(7, 7);
	add(a, 0, 0, 2);
	add(a, 0, 1, 3);

	for (i=1; i<a->nrow-1; i++) {
		add(a, i, i-1, 1);
		add(a, i, i, 2);
		add(a, i, i+1, 3);
	}

	add(a, 6, 5, 1);
	add(a, 6, 6, 2);

	printsp("a(sparse)", a);
	sptom(a, ma);
	mtxprint("a(full)", ma);
	gauss_jordan(a, b);
	printsp("solution", a);
	printv(b, a->nrow);
	deletsp(a);
}
