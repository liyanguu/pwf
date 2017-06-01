/* test_colspm.c */

#include <stdio.h>
#include <stdlib.h>
#include "colspm.h"

int main() {
	ColSpm *A;
	int i, j;

	A = colspm_init(3);
	printf("initialized\n");
	for (j = 0; j < 3; j++) {
	    	if (!colspm_addcol(A, j))
		    	printf("can't add column %d\n", j);
		printf("start adding column %d\n", j);
		colspm_add(A, j, 3.33);
		if (!colspm_add(A, j+1, 9.99))
			printf("can't added elements 9.99\n");
		printf("added elements \n");
		if (!colspm_endcol(A))
		    	printf("can't add new column to matrix\n");
		printf("finished adding column %d\n", j);
	}
	colspm_adjust(A);
	printf("finished adding elements\n");
	colspm_print(A);
	colspm_printex(A);
	colspm_clr(A);
	printf("clear the matrix\n");
	for (j = 0; j < 3; j++) {
	    	colspm_addcol(A, j);
		for (i = 0; i < 3; i++)
		    	colspm_add(A, i, 1.11);
		colspm_endcol(A);
	}
	colspm_adjust(A);
	colspm_print(A);
	colspm_printex(A);

	colspm_clr(A);
	A = colspm_del(A);
	A = colspm_init(2000);
	for (j = 0; j < 2000; j++) {
	    	if (!colspm_addcol(A, j))
		    	break;
		colspm_add(A, j, 0.99);
		colspm_endcol(A);
	}
	colspm_adjust(A);
	colspm_printex(A);

	A = colspm_del(A);
	A = colspm_init(NONZ_BUFSZ);
	for (j = 0; j < NONZ_BUFSZ; j++) {
	    	if (!colspm_addcol(A, j))
		    	break;
		colspm_add(A, j, 0.99);
		colspm_endcol(A);
	}
	colspm_adjust(A);
	colspm_printex(A);
	A = colspm_del(A);
}
