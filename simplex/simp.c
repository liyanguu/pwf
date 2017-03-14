#include <stdio.h>
#include <stdlib.h>
#include "mtx.h"
#include "simp.h"

LP makelp(Size m, Size n, Elm *cons, Elm *prof, Elm *solu) {
	LP lp;
	Size i, j;

	if ((lp = (LP) malloc(sizeof(struct lpstru))) == NULL)
		return NULL;
	lp->nrow = m;
	lp->ncol = n;
	lp->flag = UNDT;
	lp->cm = mtxalloc(m, n);
	lp->pv = mtxalloc(n, 1);
	lp->sv = mtxalloc(n, 1);
	lp->obj = mtxalloc(1, 1);
	lp->dlt_tr = mtxalloc(1, n-m);
	lp->yv = mtxalloc(m, 1);
	lp->tbm = mtxalloc(1, m);
	lp->tnbm = mtxalloc(1, n-m);
	lp->tcm = mtxalloc(m, 1);

	lp->bm_in = mtxalloc(m, m);
	for (i = 0; i < m; i++) {
		lp->bm_in->val[i][i] = 1;
		for (j = i+1; j < m; j++)
			lp->bm_in->val[i][j] = lp->bm_in->val[j][i] = 0;
	}

	lp->idx = (Size *) malloc(sizeof(Size) * n);
	for (i = 0; i < n; i++)
		lp->idx[i] = i;

	mtxread(lp->cm, cons);
	mtxread(lp->pv, prof);
	mtxread(lp->sv, solu);
	lp->bsv = rowsub(lp->sv, &lp->idx[n-m], m);
	lp->bpv = rowsub(lp->pv, &lp->idx[n-m], m);
	lp->npv = rowsub(lp->pv, lp->idx, n-m);
	lp->nbm = colsub(lp->cm, lp->idx, n-m);
	lp->bm = colsub(lp->cm, &lp->idx[n-m], m);

	return lp;
}

void printlp(LP lp) {
	Size i;

	mtxprint("constraint matrix", lp->cm);
	mtxprint("base matrix", lp->bm);
	mtxprint("none-base matrix", lp->nbm);
	mtxprint("base inverse matrix", lp->bm_in);
	xmtxprint("base solution vector", lp->bsv);
	printf("object: %10.5f\n", object(lp));
	printf("solution: \n");
	for (i = 0; i < lp->ncol - lp->nrow; i++)
		printf("%10.5f ", lp->sv->val[i][0]);
	printf("\n");
	printf("profit: \n");
	for (i = 0; i < lp->ncol - lp->nrow; i++)
		printf("%10.5f ", lp->pv->val[i][0]);
	printf("\n");
}

int simplex(LP lp) {
	Size i, j, k, r;
	Size n, m;
	Elm max, val;
	Mtx tbpv, tnpv;
	
	m = lp->nrow;
	n = lp->ncol;

	printf("multiplying bm_in... \n");
	xmtxmul(lp->tbm, lp->bpv, lp->bm_in);
	printf("multiplying nbm... \n");
	mtxmul(lp->tnbm, lp->tbm, lp->nbm);
	printf("differecing tnbm... \n");
	xmtxdif(lp->dlt_tr, lp->npv, lp->tnbm);
	xmtxmul(lp->obj, lp->bpv, lp->bsv);

	max = 0;
	for (k = j = 0; j < lp->dlt_tr->ncol; ++j)
		if (lp->dlt_tr->val[0][j] > max) {
			k = j;
			max = lp->dlt_tr->val[0][j];
		}
	if (max <= 0)
		return lp->flag = OPTM;

	printf("obtaining tcm... \n");
	colcopy(lp->tcm, lp->cm, &lp->idx[k], 1);
	mtxmul(lp->yv, lp->bm_in, lp->tcm);
	max = 0;
	for (r = i = 0; i < lp->yv->nrow; ++i) {
		val = lp->yv->val[i][0] / lp->bsv->val[i][0];
		if (val > max) {
			max = val;
			r = i;
		}
	}
	if (max == 0)
		return lp->flag = UNBD;

	/* update base inverse matrix */
	for (j = 0; j < lp->bm_in->ncol; ++j)
		lp->bm_in->val[r][j] /= lp->yv->val[r][0];
	for (i = 0; i < lp->bm_in->nrow; ++i) {
		if (i == r)
			continue;
		for (j = 0; j < lp->bm_in->ncol; ++j)
			lp->bm_in->val[i][j] -=
		lp->yv->val[i][0] * lp->bm_in->val[r][j];
	}

	/* update base solution and solution */
	swapt(Size, lp->idx[k], lp->idx[n-m+r]);
	for (i = 0; i < lp->bsv->nrow; ++i) {
		if (i != r)
			lp->bsv->val[i][0] -= lp->yv->val[i][0] / max;
		else
			lp->bsv->val[i][0] = 1 / max;
		lp->sv->val[lp->idx[n-m+i]][0] = lp->bsv->val[i][0];
	}
	lp->sv->val[lp->idx[k]][0] = 0;

	/* update others */
	rowcopy(lp->bpv, lp->pv, &lp->idx[n-m], m);
	rowcopy(lp->npv, lp->pv, lp->idx, n-m);
	colcopy(lp->nbm, lp->cm, lp->idx, n-m);
	colcopy(lp->bm, lp->cm, &lp->idx[n-m], m);

	return lp->flag = UNDT; 
}
