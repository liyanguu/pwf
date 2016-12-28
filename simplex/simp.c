#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include "mtx.h"

#define ELMMAX	DBL_MAX
#define EPS	DBL_EPSILON
#define LIM 10 /* iteration limit */
#define N 8    /* number of vars in modified problem */
#define M 3    /* number of base vars/constraints */ 
#define L 3    /* number of artificial vars */

int revsimplex(const Size n, const Size m, const Size p, const int itrlim,
	Elm *consm, Elm *prof, Elm *obj, Elm *solu);
int simplex(const Size n, const Size m, const int itrlim,
	Elm *consm, Elm *profit, Size *idx, Mtx bm_in,
	Elm *object, Elm *solu);

main() {
	Elm consm[] = {
		1, 0, 1, 0, 0, 1, 0, 0,
		1, 0, 0, -1,0, 0, 1, 0,
		0, 1, 0, 0, 1, 0, 0, 1
	};
	Elm profit[N] = {1, 0, 0, 0, 0, 0, 0, 0};
	Elm s[N] = { 0, 0, 0, 0, 0, 3, 1, 2};
	Elm maxval;
	int i, state;

	state = revsimplex(N, M, L, LIM, consm, profit, &maxval, s);
	printf("%.5f\n", (float) maxval);
	for (i = 0; i < N; ++i)
		printf(i < N-1 ? "%.5f " : "%.5f\n", s[i]);
	return state;
}

int revsimplex(const Size n, const Size m, const Size p, const int itrlim,
	Elm *consm, Elm *prof, Elm *pobj, Elm *solu) {

	int s; /* state flag */
	Mtx cm2, bm_in;
	Elm obj1, *prof1, *consm2; 
	Size i, j, *idx;

	/* phase I */
	printf("phase I\n");
	bm_in = mtxalloc(m, m);
	prof1 = malloc(sizeof(Elm) * n); /* aux profit (... -1) */ 
	idx = malloc(sizeof(Size) * n);

	for (i = 0; i < n - p; ++i)
		prof1[i] = 0;
	for ( ; i < n; ++i)
		prof1[i] = -1.0;	
	for (i = 0; i < n; ++i)
		idx[i] = i;	
	for (i = 0; i < bm_in.nrow; ++i)
		for (j = 0; j < bm_in.ncol; ++j)
			bm_in.val[i][j] = i == j ? 1.0 : 0;

	s = simplex(n, m, itrlim, consm, prof1, idx, bm_in, &obj1, solu);
	if (s < 0 || obj1 < 0) /* no feasible solution */
		return -1;

	/* phase II */
	printf("phase II\n");
	cm2 = mtxalloc(m, n - p); /* reduced constraint mtx */
	consm2 = malloc(sizeof(Elm) * m * (n-p)); /* reduced constraint */
	for (i = j = 0; i < n; ++i) /* remove artificial vars */
		if (idx[i] < n - p)	
			idx[j++] = idx[i];
	mtxread(cm2, consm, n);
	mtxwrite(cm2, consm2);	
	s = simplex(n-p, m, itrlim, consm2, prof, idx, bm_in, pobj, solu);	
	return s;
}

int simplex(const Size n, const Size m, const int itrlim,
	Elm *consm, Elm *profit, Size *idx, Mtx bm_in,
	Elm *object, Elm *solu) {
	
	int counter;
	Size i, j, k, r;
	Elm max, val;
	Mtx bm;         /* base matrix */
	Mtx cm; 	/* constraint matrix */
	Mtx pv; 	/* profit vector */
	Mtx sv; 	/* solution vector */
	Mtx nbm; 	/* none base matrix */
	Mtx bsv;  	/* base solution vector*/
	Mtx bpv, nbpv;  /* base profit vector and complement */
	Mtx dlt_tr, yv;	/* temp value */
	
	/* 1. initialization data: */
	cm = mtxalloc(m, n);
	bm = mtxalloc(m, m);
	pv = mtxalloc(n, 1);
	sv = mtxalloc(n, 1);
	nbm = mtxalloc(m, n - m);
	bsv = mtxalloc(m, 1);
	bpv = mtxalloc(m, 1);
	nbpv = mtxalloc(n - m, 1);

	mtxread(cm, consm, n);
	mtxread(pv, profit, 1);
	mtxread(sv, solu, 1);
	mtxcpy(bsv, rowsub(sv, m, &idx[n-m]));
	mtxcpy(bpv, rowsub(pv, m, &idx[n-m]));
	mtxcpy(nbpv, rowsub(pv, n - m, idx));
	mtxcpy(nbm, colsub(cm, n - m, idx));
	mtxcpy(bm, colsub(cm, m, &idx[n-m]));

	/* 2. iteration */
	for (counter = 0; counter < itrlim; ++counter) {
		mtxprint("base", bm);
		dlt_tr = mtxdif(trans(nbpv),
			mtxmul(mtxmul(trans(bpv), bm_in), nbm));
		max = 0;
		for (k = j = 0; j < dlt_tr.ncol; ++j)
			if (dlt_tr.val[0][j] > max) {
				k = j;
				max = dlt_tr.val[0][j];
			}
		if (max <= EPS) {
			mtxwrite(sv, solu);
			mtxwrite(mtxmul(trans(bpv), bsv), object);
			return 0;
		}
		yv = mtxmul(bm_in, colsub(cm, 1, &idx[k]));	
		max = 0;
		for (r = i = 0; i < yv.nrow; ++i) 	
			if ((val=yv.val[i][0]/bsv.val[i][0]) > max) {
				max = val;
				r = i;
			}
		if (max == 0) {
			printf("unbounded\n");
			return -1;
		}

		/* update */
		swap(Size, idx[k], idx[n-m+r]);
		for (j = 0; j < bm_in.ncol; ++j)
			bm_in.val[r][j] /= yv.val[r][0];
		for (i = 0; i < bm_in.nrow; ++i) {
			if (i == r)
				continue;
			for (j = 0; j < bm_in.ncol; ++j)
				bm_in.val[i][j]-=yv.val[i][0]*bm_in.val[r][j];
		}
		for (i = 0; i < bsv.nrow; ++i) {
			if (i != r)
				bsv.val[i][0] -= yv.val[i][0] / max;
			else
				bsv.val[i][0] = 1 / max;
			sv.val[idx[n-m+i]][0] = bsv.val[i][0];
		}
		sv.val[idx[k]][0] = 0;
		mtxcpy(bpv, rowsub(pv, m, &idx[n-m]));
		mtxcpy(nbpv, rowsub(pv, n - m, idx));
		mtxcpy(nbm, colsub(cm, n - m, idx));
		mtxcpy(bm, colsub(cm, m, &idx[n-m]));
		mtxprint("base solu", bsv);
	}
	printf("iteration limit reached\n");
	mtxwrite(sv, solu);
	mtxwrite(mtxmul(trans(bpv), bsv), object);
	return -1;
}

