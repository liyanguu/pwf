/* test_sp.c - basic test for spm.c */
#include <stdio.h>
#include "spm.h"
#include "mtx.h"

int main() {
	Spm a;
	int i;
	Size nz;
	Size ai[20];
	Size ap[7+1];
	Elm ax[20];

	a = spcreat(7, 7);
	spadd(a, 0, 0, 2);
	spadd(a, 0, 1, 3);

	for (i=1; i<a->nrow-1; i++) {
		spadd(a, i, i-1, 1);
		spadd(a, i, i, 2);
		spadd(a, i, i+1, 3);
	}

	spadd(a, 6, 5, 1);
	spadd(a, 6, 6, 2);

	spprint("a, (sparse)", a);
	sprem(a, 6, 6);
	sprem(a, 0, 0);
	sprem(a, 2, 2);
	sprem(a, 4, 4);
	sprem(a, 4, 0);
	spadd(a, 4, 4, 8.9);
	spprint("a, (sparse)", a);
	spcmp(a, &nz, ap, ai, ax);
	arrayprint(ap, 8);
	arrayprint(ai, 20);
	arrayprint(ax, 20);
	spdelet(a);
	return 0;
}
