#include <stdio.h> 
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>	/* for fabs() */
#include <lapacke.h>
#include <cblas.h>
#include "mtx.h"
#include "msg.h"

Size norm(Elm *v, Size lim, Elm *nr) {
	Size i, max;

	max = *nr = 0;
	for (i = 0; i < lim; i++)
		if (*nr < fabs(v[i])) {
			*nr = fabs(v[i]);
			max = i;
		}
	return max;
}

void arraycpy(Elm *dst, Elm *src, Size lim) {
	while (lim-- > 0)
		*dst++ = *src++;
}

Elm getel(Mtx m, Size i, Size j) {
	return (*m->val)[j + i * m->ncol];
}

void setel(Mtx m, Size i, Size j, Elm v) {
	(*m->val)[j + i * m->ncol] = v;
}

/* mtxread: read from element array val to matrix m */
void mtxread(Mtx m, Elm *val) {
	arraycpy(m->val[0], val, m->ncol * m->nrow);
} 

void mtxnread(Mtx m, Elm *val, int n) {
	Size i, j;

	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < m->ncol; ++j)
			setel(m, i, j, val[j + i * n]);
}

/* mtxwrite: write matrix m to element array val */
void mtxwrite(Mtx m, Elm *val) {
	Elm *p = m->val[0];
	Size nelm = m->ncol * m->nrow;
	
	arraycpy(val, p, nelm);
}

/* mtxcopy() - copy src's val to dst */
Mtx mtxcopy(Mtx dst, Mtx src) {
	Size lim = min(dst->nrow * dst->ncol, src->nrow * src->ncol);

	arraycpy(dst->val[0], src->val[0], lim);
	return dst;
}

Mtx mtxdup(Mtx a) {
	Mtx b = mtxalloc(a->nrow, a->ncol);

	if (b != NULL)
		mtxcopy(b, a);
	return b;
}


/* mtxmul -  return m = m1 * m2 */ 
Mtx mtxmul(Mtx m, Mtx m1, Mtx m2) {
	const CBLAS_ORDER order = CblasRowMajor;
	const CBLAS_TRANSPOSE trm1 = CblasNoTrans;
	const CBLAS_TRANSPOSE trm2 = CblasNoTrans;
	const Size M = m1->nrow;
	const Size N = m2->ncol;
	const Size K = m1->ncol;
	const Elm alpha = 1;
	const Elm beta = 0;
	const Size lda = m1->ncol;
	const Size ldb = m2->ncol;
	const Size ldc = m->ncol;

	cblas_dgemm(order, trm1, trm2, M, N, K, alpha, *m1->val, lda,
		*m2->val, ldb, beta, *m->val, ldc);
	return m;
}

/* xmtxmul - return m = trans(m1) * m2 */
Mtx xmtxmul(Mtx m, Mtx m1, Mtx m2) {
	const CBLAS_ORDER order = CblasRowMajor;
	const CBLAS_TRANSPOSE trm1 = CblasTrans;
	const CBLAS_TRANSPOSE trm2 = CblasNoTrans;
	const Size M = m1->nrow;
	const Size N = m2->ncol;
	const Size K = m1->ncol;
	const Elm alpha = 1;
	const Elm beta = 0;
	const Size lda = m1->ncol;
	const Size ldb = m2->ncol;
	const Size ldc = m->ncol;

	cblas_dgemm(order, trm1, trm2, M, N, K, alpha, *m1->val, lda,
		*m2->val, ldb, beta, *m->val, ldc);
	return m;
}

/* mtxdif: return m1 - m2 */
Mtx mtxdif(Mtx m, Mtx m1, Mtx m2) {
	Size i, j;
	Elm dif;

	for (i = 0; i < m1->nrow; ++i)
		for (j = 0; j < m1->ncol; ++j) {
			dif = getel(m1, i, j) - getel(m2, i, j);
			setel(m, i, j, dif);
		}
	return m;
}

Mtx xmtxdif(Mtx m, Mtx m1, Mtx m2) {
	Size i, j;
	Elm dif;

	for (i = 0; i < m2->nrow; ++i)
		for (j = 0; j < m2->ncol; ++j) {
			dif = getel(m1, j, i) - getel(m2, i, j);
			setel(m, i, j, dif);
		}
	return m;
}

Mtx trans(Mtx m) {
	Size i, j;
	Elm e;
	Mtx transm = mtxalloc(m->ncol, m->nrow);

	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < m->ncol; ++j) {
			e = getel(m, i, j);
			setel(transm, i, j, e);
		}
	return transm;
}

/* return submatrix ms = m(row index, :) */
Mtx rowsub(Mtx m, Size *idx, Size nrow) {
	Size i, j;
	Mtx ms;
	Elm e;

	ms = mtxalloc(nrow, m->ncol);
	for (i = 0; i < nrow; ++i)
		for (j = 0; j < m->ncol; ++j) {
			e = getel(m, idx[i], j);
			setel(ms, i, j, e);
		}
	return ms;
}

/* copy row sub-matrix of m to mc, return mc */
Mtx rowcopy(Mtx mc, Mtx m, Size *idx, Size nrow) {
	Size i, j;
	Elm e;

	for (i = 0; i < nrow; ++i)
		for (j = 0; j < m->ncol; ++j) {
			e = getel(m, idx[i], j);
			setel(mc, i, j, e);
		}
	return mc;
}

/* submatrix m(:, col index) */
Mtx colsub(Mtx m, Size *idx, Size ncol) {
	Size i, j;
	Mtx ms;
	Elm e;

	ms = mtxalloc(m->nrow, ncol);
	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < ncol; ++j) {
			e = getel(m, i, idx[j]);
			setel(ms, i, j, e);
		}
	return ms;
}

/* copy column sub-matrix of m to mc, return mc */
Mtx colcopy(Mtx mc, Mtx m, Size *idx, Size ncol) {
	Size i, j;
	Elm e;

	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < ncol; ++j) {
			e = getel(m, i, idx[j]);
			setel(mc, i, j, e);
		}
	return mc;
}

/* mtxprint() - print a matrix */
void mtxprint(char *ttl, Mtx a) {
	Size i, j;
	
	printf("%s\n", ttl);
	for (i = 0; i < a->nrow; ++i) {
		for (j = 0; j < a->ncol; ++j)
			printf("%10.5f ", getel(a, i, j)); 
		putchar('\n');
	}
	putchar('\n');
}

void xmtxprint(char *ttl, Mtx a) {
	Size i, j;

	printf("%s\n", ttl);
	for (i = 0; i < a->ncol; ++i) {
		for (j = 0; j < a->nrow; ++j)
			printf("%10.5f ", getel(a, j, i));
		putchar('\n');
	}
	putchar('\n');
}

void mtxnprt(char *ttl, ...) {
	va_list ap;
	char *p;

	va_start(ap, ttl);
	for (p = ttl; *p; p++) {
		if (*p != '%') {
			putchar(*p);
			continue;
		}
		mtxprint("", va_arg(ap, Mtx));
	}
	va_end(ap);
} 

Mtx reshape(Mtx m, Size newr, Size newc) {
	Mtx newm;
	Size newi, newj, i, j, k, si, sj;
	Elm e;

	if (newr == m->nrow && newc == m->ncol)
		return m;
	newm = mtxalloc(newr, newc);
	if (m->nrow >= newr) {
		i = m->nrow - newr;
		newi = 0;
	} else {
		i = 0;
		newi = newr - m->nrow;
		for (si = 0; si < newi; si++)
			for (sj = 0; sj < newc; sj++)
				setel(newm, si, sj, 0);
	}
	if (m->ncol >= newc) {
		j = m->ncol - newc;
		newj = 0;
	} else {
		j = 0;
		newj = newc - m->ncol;
		for (si = newi; si < newr; si++)
			for (sj = 0; sj < newj; sj++)
				setel(newm, si, sj, 0);
	}
	for (si = newi; si < newr; si++, i++)
		for (sj = newj, k = j; sj < newc; sj++, k++) {
			e = getel(m, i, k);
			setel(newm, si, sj, e);
		}
	mtxfree(m);
	return newm;
}

/* mtxalloc - allocate and return a nrow by ncol matrix.
	if there is not enough memory, return NULL */
Mtx mtxalloc(Size nrow, Size ncol) {
	Size i;
	Mtx tmp;
	Elm *p;

	if ((tmp = (Mtx) malloc(sizeof(*tmp))) == NULL 
	|| (tmp->val = (Rv) malloc(sizeof(Elm *) * nrow)) == NULL
	|| (p = (Elm *) malloc(sizeof(Elm) * nrow * ncol)) == NULL) {
		msg(stderr, "mtx: mtxalloc: not enough memory\n");
		return NULL;
	}
	tmp->nrow = nrow;
	tmp->ncol = ncol;
	for (i = 0; i < nrow; i++)
		tmp->val[i] = p + i * ncol;
	return tmp;
}

/* mtxfree() - free mtx m */ 
void mtxfree(Mtx m) {
	if (m == NULL)
		return;
	free(*m->val);
	free(m->val);
	free(m);
}

/* exrow - exchange row p and q of matrix m */
Mtx exrow(Mtx m, Size p, Size q) {
	Size i;
	Elm tmp;

	for (i = 0; i < m->ncol; ++i) {
		tmp = getel(m, p, i);
		setel(m, p, i, getel(m, q, i));
		setel(m, q, i, tmp);
	}
	return m;
}

Mtx excol(Mtx m, Size p, Size q) {
	Elm tmp;
	Size i;

	for (i = 0; i < m->nrow; ++i) {
		tmp = getel(m, i, p);
		setel(m, i, p, getel(m, i, q));
		setel(m, i, q, tmp);
	}
	return m;
}

/* read a column to matrix */
void readcol(Mtx m, Size colno, Elm *colv) {
	Size i;

	for (i = 0; i < m->nrow; ++i) {
		setel(m, i, colno, colv[i]);
	}
}

/* solve linear equations AX = B */
int mtxsv(Mtx A, Mtx b) {
	Size *ipiv = calloc(A->nrow, sizeof(Size));
	Size n = A->nrow;
	Size nrhs = b->ncol;
	Size lda = A->ncol;
	Size ldb = b->ncol;
	int info;

	info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs,
		*A->val, lda, ipiv, *b->val, ldb);
	free(ipiv);
	return info;
}

Mtx mtxinv(Mtx m) {
	Mtx mcopy = mtxdup(m);
	Size N = m->nrow;
	Size lda = m->ncol;
	Size *ipiv = malloc(m->nrow * sizeof(Size));
	int flag;

	flag = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N,
		*mcopy->val, lda, ipiv);
	if (flag) {
		free(ipiv);
		return NULL;
	}
	flag = LAPACKE_dgetri(LAPACK_ROW_MAJOR, N,
		*mcopy->val, lda, ipiv);
	free(ipiv);
	if (flag)
		return NULL;
	else
		return mcopy;
}

Mtx mtxinv2(Mtx m) {
	Mtx mcopy = mtxdup(m);
	Mtx minv = mtxalloc(m->nrow, m->ncol);
	Size i, j;
	int flag;

	for (i = 0; i < m->nrow; i++)
		for (j = 0; j < m->ncol; j++)
			setel(m, i, j, (i == j) ? 1 : 0);
	flag = mtxsv(mcopy, minv);
	mtxfree(mcopy);
	if (flag == 0)
		return minv;
	else
		return NULL;
}

int ludcmp(Mtx a, Size *indx, Elm *d) {
	Elm *vv = (Elm *) malloc(sizeof(Elm) * a->nrow);	
	Size i, j, k, imax;
	Elm aamax, sum, dum;

	for (i = 0; i < a->nrow; ++i) {
		aamax = 0;
		for (j = 0; j < a->ncol; ++j)
			if (aamax < fabs(a->val[i][j]))
				aamax = fabs(a->val[i][j]);
		if (aamax == 0) {
			msg(stderr, "mtx: ludcmp: singular matrix\n");
			return 0;
		}
		vv[i] = 1. /aamax;
		msg(stderr, "mtx: ludcmp: vv[%4d]=%10.4f\n", i, vv[i]);
	}
	for (j = 0; j < a->ncol; ++j) {
		for (i = 0; i < j; ++i) {
			sum = a->val[i][j];
			for (k = 0; k < i; ++k)
				sum -= a->val[i][k]*a->val[k][j];
			a->val[i][j] = sum;
		}
		imax = aamax = 0;
		for (i = j; i < a->nrow; ++i) {
			sum = a->val[i][j];
			for (k = 0; k < j; ++k)
				sum -= a->val[i][k]*a->val[k][j];
			a->val[i][j] = sum;
			dum = vv[i] * fabs(sum);
			if (dum >= aamax) {
				imax = i;
				aamax = dum;
			}
		}
		if (j != imax) {
			exrow(a, j, imax);
			*d *= -1;	
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a->val[j][j] == 0)
			a->val[i][j] = ZEROEPS;
		dum = 1. / a->val[j][j];
		for (i = j+1; i < a->nrow; ++i)
			a->val[i][j] *= dum;
	}
	free(vv);
	return 1;
}

void lubksb(Mtx a, Size *indx, Elm *b) {
	Size ii, ll, i, j;
	Elm sum;

	ii = -1;
	for (i = 0; i < a->nrow; ++i) {
		ll = indx[i];
		sum = b[ll];
		b[ll] = b[i];
		if (ii != -1)
			for (j = ii; j < i; ++j)
				sum -= a->val[i][j]*b[j];
		else if (sum != 0)
			ii = i;
		b[i] = sum;
	}
	for (i = a->nrow-1; i >= 0; --i) {
		sum = b[i];
		for (j = i+1; j < a->ncol; ++j)
			sum -= a->val[i][j]*b[j];
		b[i] = sum/a->val[i][i];
	}
}

/* return the inverse matrix of a, a unchanged */
Mtx inverse(Mtx a) {
	Mtx trinva, inv, acopy;
	Elm det;
	Size *indx = (Size *) malloc(sizeof(Size) * a->nrow);
	Size i, j;

	acopy = mtxdup(a);
	trinva = mtxalloc(a->nrow, a->ncol);
	ludcmp(acopy, indx, &det);
	for (i = 0; i < trinva->nrow; i++) {
		for (j = 0; j < trinva->ncol; j++)
			trinva->val[i][j] = .0;
		trinva->val[i][i] = 1.0;
		lubksb(acopy, indx, trinva->val[i]);
	}
	inv = trans(trinva);
	free(indx);
	mtxfree(trinva);
	printf("freed trinva\n");
	mtxfree(acopy);
	printf("freed acopy\n");
	return inv;
}

#define EPS 1E-6
#define NMAX 500
enum { RESTART, SINGULAR, DONE, OVERLIMIT };

Elm g[NMAX];
Elm h[NMAX];
Elm xi[NMAX];
Elm xj[NMAX];
Elm anum, aden, rp, bsq, gg, dgg, gam, eps2;

int sparsecalc(Elm *b, Size n, void (*asub)(Elm *, Elm *),
		void (*atsub)(Elm*, Elm*), Elm *x, Elm *rsq) {
	int j, iter;

	(*asub)(x, xi);
	rp = .0;
	bsq = .0;
	for (j = 0; j < n; j++) {
		bsq += b[j] * b[j];
		xi[j] = xi[j] - b[j];
		rp += xi[j] * xi[j];
	}
	(*atsub)(xi, g);
	for (j = 0; j < n; j++) {
		g[j] = -g[j];
		h[j] = g[j];
	}
	for (iter = 1; iter <= 10 * n; iter++) {
		(*asub)(h, xi);
		anum = .0;
		aden = .0;
		for (j = 0; j < n; j++) {
			anum += g[j]*h[j];
			aden += xi[j] * xi[j];
		}
		if (aden == .0)
			return SINGULAR;
		anum /= aden;
		for (j = 0; j < n; j++) {
			xi[j] = x[j];
			x[j] += anum * h[j];
		}
		(*asub)(x, xj);
		*rsq = .0;
		for (j = 0; j < n; j++) {
			xj[j] -= b[j];
			*rsq += xj[j] * xj[j];
		}
		if (*rsq == rp || *rsq <= bsq*eps2)
			return DONE;
		if (*rsq > rp) {
			for (j = 0; j < n; j++) 
			x[j] = xi[j];
			return RESTART;
		}
		rp = *rsq;
		(*atsub)(xj, xi);
		gg = .0;
		dgg = .0;
		for (j = 0; j < n; j++) {
			gg += g[j] * g[j];
			dgg += (xi[j] + g[j]) * xi[j];
		}
		if (gg == .0)
			return DONE;
		gam = dgg/gg;
		for (j = 0; j < n; j++) {
			g[j] = -xi[j];
			h[j] = g[j] + gam * h[j];
		}
	}
	printf("too many iterations\n");
	return OVERLIMIT;
}

void sparse(Elm *b, Size n, void (*asub)(Elm *, Elm *),
	void (*atsub)(Elm*, Elm*), Elm *x, Elm *rsq) {
	Size irst;
	eps2 = n * EPS * EPS;

	for (irst = 1; irst < 3; irst++)
		switch (sparsecalc(b, n, asub, atsub, x, rsq)) {
		case DONE:
			return;
		case RESTART:
			continue;
		case SINGULAR:
			msg(stderr, "mtx: sparse: very singular matrix\n");
			return;
		case OVERLIMIT:
			msg(stderr, "mtx: sparse: iteration limit reached\n");
			return;
		}
}
