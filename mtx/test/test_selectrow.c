/* test_selectrow.c - select */

#include <stdio.h>
#include "mtx.h"

void printarray(Elm *v, int len);

int main() {
	Mtx a = mtxalloc(10, 10);
	Elm buf[5];
	Size index[5] = { 1, 3, 5, 7, 9 };
	int i, j;

	for (i = 0; i < a->nrow; i++)
		for (j = 0; j < a->nrow; j++)
			a->val[i][j] = i + j;
	mtxprint("A", a);
	selectrow(a, 0, buf, 5, index); 
	printarray(buf, 5);
	selectrow(a, 3, buf, 5, index); 
	printarray(buf, 5);
	return 0;
}

void printarray(Elm *v, int len) {
	while (len-- > 0)
		printf(" %10.4f", *v++);
	printf("\n");
}
