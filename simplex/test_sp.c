#include <stdio.h>
#include "spm.h"

int main() {
	Spm a;
	int i, j;
	Elm b[] = { 5, 6, 6, 6, 6, 6, 3 };

	a = spcreat(7, 7);
	spinsert(a, 0, 0, 2);
	spinsert(a, 0, 1, 3);

	for (i=1; i<a->nrow-1; i++) {
		spinsert(a, i, i-1, 1);
		spinsert(a, i, i, 2);
		spinsert(a, i, i+1, 3);
	}

	spinsert(a, 6, 5, 1);
	spinsert(a, 6, 6, 2);

	spprint("a, (sparse)", a);
	printv(b, a->nrow);
	gauss_jordan(a, b);
	spprint("a, reduced to diagnal", a);
	printv(b, a->nrow);
	spdelet(a);
}
