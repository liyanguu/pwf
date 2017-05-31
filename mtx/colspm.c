/* colspm.c
 * 描述：压缩主列式(compressed column major)稀疏矩阵的实现
 * 2017-5-30
 */
#include <stdio.h> 
#include <stdlib.h>
#include "colspm.h"

ColSpm *colspm_alloc(void) {
    	return malloc(sizeof (struct colspm));
}

ColSpm *colspm_init(Size dim, Size nonz) {
    	ColSpm *matx;

	if ((matx = colspm_alloc()) == NULL)
	    	return NULL;
	matx->elms = malloc(nonz * sizeof (Elm));
	matx->pos = malloc(dim * sizeof (Size) + 1);
	matx->rows = malloc(nonz * sizeof (Size));
	matx->colbuf = malloc(dim * sizeof (Elm));
	if (matx->colbuf == NULL)
	    	return NULL;
	matx->pos[0] = 0;
	matx->nonz = nonz;
	matx->dim = dim;
	matx->cnt = 0;
	return matx;
}

ColSpm *colspm_del(ColSpm *matx) {
    	if (matx == NULL)
	    	return NULL;
	free(matx->elms);
	free(matx->pos);
	free(matx->rows);
	free(matx->colbuf);
	free(matx);
	return NULL;
}

ColSpm *colspm_clr(ColSpm *matx) {
    	matx->cnt = 0;
	return matx;
}

int colspm_addcol(ColSpm *matx, Size colno) {
    	if (colno >= matx->dim)
	    	return 0;
    	matx->colcnt = 0;
	matx->curcol = colno;
	return 1;
}

int colspm_add(ColSpm *matx, Size rowno, Elm elm) {
    	Size cnt = matx->cnt;

	if (cnt >= matx->nonz || matx->colcnt >= matx->dim 
	|| rowno >= matx->dim)
	    	return 0;
    	matx->elms[cnt] = elm;
    	matx->rows[cnt] = rowno;
	matx->cnt++;
	matx->colcnt++;
	return 1;
}

void colspm_endcol(ColSpm *matx) {
    	matx->pos[matx->curcol + 1] = matx->cnt;
}

void colspm_print(ColSpm *matx) {
    	Size i, colno, rowno, rowlim;

	printf("%d, %d\n", matx->dim, matx->nonz);
	for (colno = 0; colno < matx->dim; colno++) {
		rowlim = matx->pos[colno+1];
		rowno = matx->pos[colno];
		for (i = 0; rowno < rowlim; rowno++, i++) {
		    	for ( ; i < matx->rows[rowno]; i++)
			    	matx->colbuf[i] = 0;
			matx->colbuf[i] = matx->elms[rowno];
		}
		for ( ; i < matx->dim; i++)
			  matx->colbuf[i] = 0;
		arrayprint(matx->colbuf, matx->dim);
	}
}

void colspm_printex(ColSpm *matx) {
    	Size rowno, colno, rowlim;

	printf("%d, %d\n", matx->dim, matx->nonz);
	for (colno = 0; colno < matx->dim; colno++) {
	    	rowno = matx->pos[colno];
		rowlim = matx->pos[colno+1];
		for ( ; rowno < rowlim; rowno++)
		    	printf("%d %d %.4f\n", colno, matx->rows[rowno],
				matx->elms[rowno]);
	}
}
