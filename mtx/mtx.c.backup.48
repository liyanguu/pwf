#include <stdio.h> 
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>	/* for fabs() */
#include "mtx.h"

/* mtxalloc - allocate and return a nrow by ncol matrix.
	if there is not enough memory, return NULL */
Mtx mtxalloc(Size nrow, Size ncol) {
	Size i;
	Elm *p;
	Mtx tmp;

	if (nrow <= 0 || nrow <= 0
		|| (tmp = (Mtx) malloc(sizeof(struct mtx))) == NULL
		|| (tmp->val = (Rv) malloc(sizeof(Elm*) * nrow)) == NULL
		|| (p = (Elm *) malloc(sizeof(Elm) * nrow * ncol)) == NULL)
		return NULL;
	tmp->nrow = nrow;
	tmp->ncol = ncol;
	for (i = 0; i < nrow; i++)
		tmp->val[i] = p + ncol * i;
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

Size norm(Elm *v, Size lim, Elm *nr) {
	Size i, max;

	*nr = 0;
	for (i = 0; i < lim; i++)
		if (*nr < fabs(v[i])) {
			*nr = fabs(v[i]);
			max = i;
		}
	return i;
}

/* arraycpy - copy src to dst */
void arraycpy(Elm *dst, Elm *src, Size lim) {
	while (lim-- > 0)
		*dst++ = *src++;
}

/* mtxread - read from element array val to matrix m */
void mtxread(Mtx m, Elm *val) {
	arraycpy(*m->val, val, m->ncol * m->nrow);
} 

/* mtxnread - seem as mtxread, but the leading dimision of val is n */
void mtxnread(Mtx m, Elm *val, int n) {
	Size i, j;

	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < m->ncol; ++j)
			m->val[i][j] = val[j + i * n];
}

/* mtxwrite - write matrix m to array val */
void mtxwrite(Mtx m, Elm *val) {
	arraycpy(val, *m->val, m->ncol * m->nrow);
}

/* mtxcopy() - copy m2's val to m1, m1 must have the same size of m2 */
Mtx mtxcopy(Mtx m1, Mtx m2) {
	Size i, j;

	for (i = 0; i < m1->nrow; ++i)
		for (j = 0; j < m1->ncol; ++j)
			m1->val[i][j] = m2->val[i][j]; /* a pointer */
	return m1;
}

/* mtxdup - allocate a mtx which is copy of a */
Mtx mtxdup(Mtx a) {
	Mtx b = mtxalloc(a->nrow, a->ncol);

	if (b != NULL)
		mtxcopy(b, a);
	return b;
}


/* mtxmul -  return m = m1 * m2 */ 
Mtx mtxmul(Mtx m, Mtx m1, Mtx m2) {
	Size i, j, k;

	for (i = 0; i < m1->nrow; ++i)
		for (j = 0; j < m2->ncol; ++j) {
			m->val[i][j] = 0;
			for (k = 0; k < m1->ncol; ++k)
				m->val[i][j]+=m1->val[i][k]*m2->val[k][j];
		}
	return m;
}

Mtx xmtxmul(Mtx m, Mtx m1, Mtx m2) {
	Size i, j, k;

	for (i = 0; i < m1->ncol; ++i)
		for (j = 0; j < m2->ncol; ++j) {
			m->val[i][j] = 0;
			for (k = 0; k < m1->nrow; ++k)
				m->val[i][j]+=m1->val[k][i]*m2->val[k][j];
		}
	return m;
}

/* mtxdif: return m1 - m2 */
Mtx mtxdif(Mtx m, Mtx m1, Mtx m2) {
	Size i, j;

	for (i = 0; i < m1->nrow; ++i)
		for (j = 0; j < m1->ncol; ++j)
			m->val[i][j] = m1->val[i][j] - m2->val[i][j];
	return m;
}

Mtx xmtxdif(Mtx m, Mtx m1, Mtx m2) {
	Size i, j;

	for (i = 0; i < m2->nrow; ++i)
		for (j = 0; j < m2->ncol; ++j)
			m->val[i][j] = m1->val[j][i] - m2->val[i][j];
	return m;
}

Mtx trans(Mtx m) {
	Size i, j;
	Mtx mt = mtxalloc(m->ncol, m->nrow);

	for (i = 0; i < mt->nrow; ++i)
		for (j = 0; j < mt->ncol; ++j)
			mt->val[i][j] = m->val[j][i];
	return mt;
}

/* return submatrix ms = m(row index, :) */
Mtx rowsub(Mtx m, Size *idx, Size nrow) {
	Size i, j;
	Mtx ms;

	ms = mtxalloc(nrow, m->ncol);
	for (i = 0; i < nrow; ++i)
		for (j = 0; j < m->ncol; ++j)
			ms->val[i][j] = m->val[idx[i]][j];
	return ms;
}

Mtx rowcopy(Mtx mc, Mtx m, Size *idx, Size nrow) {
	Size i, j;

	for (i = 0; i < nrow; ++i)
		for (j = 0; j < m->ncol; ++j)
			mc->val[i][j] = m->val[idx[i]][j];
	return mc;
}

/* submatrix m(:, col index) */
Mtx colsub(Mtx m, Size *idx, Size ncol) {
	Size i, j;
	Mtx ms;

	ms = mtxalloc(m->nrow, ncol);
	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < ncol; ++j)
			ms->val[i][j] = m->val[i][idx[j]];
	return ms;
}

Mtx colcopy(Mtx mc, Mtx m, Size *idx, Size ncol) {
	Size i, j;

	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < ncol; ++j)
			mc->val[i][j] = m->val[i][idx[j]];
	return mc;
}

int selectrow(Mtx m, Size rowno, Elm *val, Size len, Size *idx) {
	Size i, n;

	if (rowno < 0 || rowno >= m->nrow)
		return 0;
	for (i = 0; i < len; i++) {
		if ((n = idx[i]) < 0 || n >= m->ncol)
			return 0;
		val[i] = m->val[rowno][n];
	}
	return 1;
}

Mtx reshape(Mtx m, Size newr, Size newc) {
	Mtx newm;
	Size newi, newj, i, j, k, si, sj;

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
				newm->val[si][sj] = 0;
	}
	if (m->ncol >= newc) {
		j = m->ncol - newc;
		newj = 0;
	} else {
		j = 0;
		newj = newc - m->ncol;
		for (si = newi; si < newr; si++)
			for (sj = 0; sj < newj; sj++)
				newm->val[si][sj] = 0;
	}
	for (si = newi; si < newr; si++, i++)
		for (sj = newj, k = j; sj < newc; sj++, k++)
			newm->val[si][sj] = m->val[i][k];	
	mtxfree(m);
	return newm;
}


/* fmtxprt() - print a matrix to file fp */
void fmtxprt(FILE *fp, char *title, Mtx a) {
	Size i, j;
	
	fprintf(fp, "%s =\n", title);
	for (i = 0; i < a->nrow; ++i) {
		for (j = 0; j < a->ncol; ++j)
			fprintf(fp, "%10.5f ", a->val[i][j]);
		putc('\n', fp);
	}
	putc('\n', fp);
}

/* fxmtxprt - same as fmtxprt, but a transposed */
void fxmtxprt(FILE *fp, char *title, Mtx a) {
	Size i, j;

	fprintf(fp, "%s =\n", title);
	for (i = 0; i < a->ncol; ++i) {
		for (j = 0; j < a->nrow; ++j)
			fprintf(fp, "%10.5f ", a->val[j][i]);
		putc('\n', fp);
	}
	putc('\n', fp);
}

/* fmtxnprt - print matries according to fmt 
	a fmt is a string of :
	normal chars	-	print as is
	%[m,t]			-	matrix specifier, if m, print a matrix.
						if t, print the transpose */
void fmtxnprt(FILE *fp, char *fmt, ...) {
	va_list ap;

	va_start(ap, fmt);
	vfmtxnprt(fp, fmt, ap);
	va_end(ap);
}

void mtxnprt(char *fmt, ...) {
	va_list ap;

	va_start(ap, fmt);
	vfmtxnprt(stdout, fmt, ap);
	va_end(ap);
}

void vfmtxnprt(FILE *fp, char *fmt, va_list arg) {
	char *p;

	for (p = fmt; *p; p++) {
		if (*p != '%') {
			putc(*p, fp);
			continue;
		}
		if (*++p == 't')
			fxmtxprt(fp, "", va_arg(arg, Mtx));
		else if (*p == 'm')
			fmtxprt(fp, "", va_arg(arg, Mtx));
		else
			--p;
	}
}

/* exrow - exchange row p and q of matrix m */
Mtx exrow(Mtx m, Size p, Size q) {
	swapt(Elm *, m->val[p], m->val[q]);
	return m;
}

Mtx excol(Mtx m, Size i, Size j) {
	Size p;

	for (p = 0; p < m->nrow; ++p)
		swapt(Elm, m->val[p][i], m->val[p][j]);
	return m;
}

void readcol(Mtx m, Size colno, Elm *colv) {
	Size i;

	for (i = 0; i < m->nrow; ++i)
		m->val[i][colno] = colv[i];
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
			fprintf(stderr, "mtx: ludcmp: singular matrix\n");
			return 0;
		}
		vv[i] = 1. /aamax;
		fprintf(stderr, "mtx: ludcmp: vv[%4d]=%10.4f\n", i, vv[i]);
	}
	for (j = 0; j < a->ncol; ++j) {
		for (i = 0; i < j; ++i) {
			sum = a->val[i][j];
			for (k = 0; k < i; ++k)
				sum -= a->val[i][k]*a->val[k][j];
			a->val[i][j] = sum;
		}
		aamax = 0;
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

int chodcm(Mtx a, Elm *d, Elm *t) {
	Size i, j, k;
	Elm sum;

	for (i = 0; i < a->nrow; ++i) {
		sum = a->val[i][i];
		for (j = 0; j < i; ++j) {
			t[j] = a->val[j][i];
			for (k = 0; k < j; ++k)
				t[j] -= t[k] * a->val[j][k];
			if (d[j] == 0) {
				if (t[j] != 0) {
					fprintf(stderr, "No Cholesky Decompsition\n");
					return 0;
				} else
					a->val[i][j] = 1.0;
			} else
				a->val[i][j] = t[j] / d[j];
			sum -= t[j] * a->val[i][j];
		}
		d[j] = sum;
	}
	return 1;
}

int chobsb(Mtx a, Elm *d, Elm *b) {
	Size i, j;
	Elm sum;

	for (i = 0; i < a->nrow; ++i) {
		sum = b[i];
		for (j = 0; j < i; ++j)
			sum -= a->val[i][j] * b[j];
		b[i] = sum;
	}
	for (i = a->nrow-1; i >= 0; --i) {
		if (d[i] == 0) {
			fprintf(stderr, "singular matrix\n");
			return 0;
		}
		sum = b[i] / d[i];
		for (j = i+1; j < a->ncol; ++j)
			sum -= a->val[j][i] * b[j];
		b[i] = sum;
	}
	return 1;
}

Mtx inverse(Mtx a) {
	Mtx b, inva, inv;
	Elm det;
	Size *indx = (Size *) malloc(sizeof(Size) * a->nrow);
	Size i, j;

	b = mtxdup(a);
	ludcmp(b, indx, &det);
	inva = mtxalloc(a->nrow, a->ncol);
	for (i = 0; i < inva->nrow; i++) {
		for (j = 0; j < inva->ncol; j++)
			inva->val[i][j] = .0;
		inva->val[i][i] = 1.0;
		lubksb(b, indx, inva->val[i]);
	}
	free(indx);
	mtxfree(b);
	inv = trans(inva);
	mtxfree(inva);
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
			fprintf(stderr, "mtx: sparse: very singular matrix\n");
			return;
		case OVERLIMIT:
			fprintf(stderr, "mtx: sparse: iteration limit reached\n");
			return;
		}
}
