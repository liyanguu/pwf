/* _530.c - test colspm */

#include <stdio.h>
#include <stdlib.h>
#include "colspm.h"

int main() {
	ColSpm *A;
	int i, j;

	A = colspm_init(3, 5);
	for (j = 0; j < 3; j++) {
	    	colspm_addcol(A, j);
		colspm_add(A, j, 3.33);
		colspm_add(A, j+1, 9.99);
		colspm_endcol(A);
	}
	colspm_print(A);
	colspm_printex(A);
	colspm_clr(A);
	for (j = 0; j < 3; j++) {
	    	colspm_addcol(A, j);
		colspm_add(A, j, 1.245);
		colspm_add(A, j+1, 2.89);
		colspm_endcol(A);
	}
	colspm_print(A);
	colspm_printex(A);

	A = colspm_del(A);
}
