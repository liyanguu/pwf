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

ColSpm *colspm_init(Size dim) {
    	ColSpm *matx;

	if ((matx = colspm_alloc()) == NULL)
	    	return NULL;
	matx->elms = malloc(NONZ_BUFSZ * sizeof (Elm));
	matx->rows = malloc(NONZ_BUFSZ * sizeof (Size));
	matx->pos = malloc(dim * sizeof (Size) + 1);
	matx->index = malloc(dim * sizeof (Size));
	matx->colbuf = malloc(dim * sizeof (Elm));
	if (matx->colbuf == NULL)
	    	return NULL;
	matx->pos[0] = 0;
	matx->nz = NONZ_BUFSZ;
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
	free(matx->index);
	free(matx);
	return NULL;
}

ColSpm *colspm_clr(ColSpm *matx) {
    	matx->cnt = 0;
	return matx;
}

int colspm_addcol(ColSpm *matx, Size colno) {
    	if (colno < 0 || colno >= matx->dim)
	    	return 0;
    	matx->colcnt = 0;
	matx->curcol = colno;
	return 1;
}

int colspm_add(ColSpm *matx, Size rowno, Elm elm) {
    	Size colcnt = matx->colcnt;

	if (rowno < 0 || rowno >= matx->dim 
	|| colcnt >= matx->dim)
	    	return 0;
    	matx->colbuf[colcnt] = elm;
    	matx->index[colcnt] = rowno;
	matx->colcnt++;
	return 1;
}

int colspm_endcol(ColSpm *matx) {
    	Size i;
	Size cnt = matx->cnt;
	Elm *newelms;	
	Size *newrows; 

  	if ((cnt += matx->colcnt) > matx->nz) {
	    	matx->nz += max(matx->colcnt, NONZ_BUFSZ);
	    	if ((newelms = realloc(matx->elms, matx->nz * sizeof (Elm))) == NULL
		|| (newrows = realloc(matx->rows, matx->nz * sizeof (Size))) == NULL)
			return 0;
		matx->elms = newelms;
		matx->rows = newrows;
	}
	newelms = matx->elms + matx->cnt;
	newrows = matx->rows + matx->cnt;
	for (i = 0; i < matx->colcnt; i++) {
		newelms[i] = matx->colbuf[i];
		newrows[i] = matx->index[i];
/* 		printf("adding %d, %f\n", newrows[i], newelms[i]); */
	}
	matx->cnt = cnt;
    	matx->pos[matx->curcol + 1] = cnt;
/* 	printf("total elms: %d\n", cnt); */
	return 1;
}

void colspm_adjust(ColSpm *matx) {
/*	printf("实际元素数：%d\n", matx->cnt); 
	printf("已分配元素数：%d\n", matx->nz); */
    	if (matx->cnt == matx->nz)
	    	return;
	matx->elms = realloc(matx->elms, matx->cnt * sizeof (Elm));
	matx->rows = realloc(matx->rows, matx->cnt * sizeof (Size));
	matx->nz = matx->cnt;
}

void colspm_print(ColSpm *matx) {
    	Size i, colno, rowno, rowlim;

	printf("%d, %d\n", matx->dim, matx->nz);
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

	printf("矩阵维度：%d\n", matx->dim);
	printf("非零元数：%d\n", matx->nz);
	printf("%s %s %s\n", "列  ", "行  ", "  值");
	for (colno = 0; colno < matx->dim; colno++) {
	    	rowno = matx->pos[colno];
		rowlim = matx->pos[colno+1];
		for ( ; rowno < rowlim; rowno++)
		    	printf("%4d %4d %.4f\n", colno, matx->rows[rowno],
				matx->elms[rowno]);
	}
}
