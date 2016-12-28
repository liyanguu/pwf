#include <stdio.h> 
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>	/* for fabs() */
#include "mtx.h"
#include "../msg.h"

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

Spm spalloc(void) {
	Spm m;

	return (Spm) malloc(sizeof(*m));
}

struct spel *spelalloc(Size n) {
	return (struct spel *) malloc(sizeof(struct spel) * n);
}

Spm makesp(Size nrow, Size ncol) {
	Spm m;
	Size i;

	m = spalloc();
	m->val = spelalloc(nrow);
	for (i=0; i < nrow; i++) {
		m->val[i].v = .0;
		m->val[i].colno = i;
		m->val[i].nextcol = NULL;
	}
	m->nrow = nrow;
	m->ncol = ncol;
	return m;
}

struct spel *getpos(Spm a, Size rowno, Size pos) { 
	struct spel *p;

	if (rowno < 0 || rowno >= a->nrow) {
		msg(stderr, "mtx: getpos: out of range %d\n", rowno);
		exit(1);
	}
	for (p=&a->val[rowno]; pos-- > 0 && p != NULL; p=p->nextcol)
		;
	return p;
}

struct spel *find(Spm a, Size rowno, Size colno) {
	struct spel *p;
	int i;

	for (i=0; (p=getpos(a, rowno, i)) != NULL; i++)
		if (p->colno == colno)
			return p;
	return NULL;
}

void add(Spm a, Size rowno, Size colno, Elm v) {
	struct spel *p;

	if ((p = find(a, rowno, colno)) != NULL) {
		p->v = v;
		return;
	} else {
		p = spelalloc(1);
		p->colno = colno;
		p->v = v;
		p->nextcol = a->val[rowno].nextcol;
		a->val[rowno].nextcol = p;
		return;
	}
}

void delete(Spm a, Size rowno, Size colno) {
	struct spel *p, *q;

	p = getpos(a, rowno, 0);
	for ( ; p->nextcol != NULL; p = p->nextcol)
		if (p->nextcol->colno == colno) {
			q = p->nextcol;
			p->nextcol = q->nextcol;
			msg(stderr, "mtx: delete: deleting %d, %d\n", rowno,q->colno);
			free(q);
			return;
		}
}

void freerow(struct spel *p) {
	struct spel *q;

	for ( ; p != NULL; p = q) {
		msg(stderr, "mtx: freerow: deleting %d\n", p->colno); 
		q = p->nextcol;
		free(p);
	}
}

void deletsp(Spm a) {
	struct spel *p;
	Size i, j, nrow;

	for (i = 0; i < a->nrow; i++)
		freerow(a->val[i].nextcol);
	free(a->val);
	free(a);
}


Elm buf[100];

void printsp(char *title, Spm a) {
	struct spel *p;
	Size i, j;

	printf("%s =\n", title);
	for (i=0; i < a->nrow; i++) {
		for (j = 0; j < a->ncol; j++)
			buf[j] = .0;
		for (j = 0; (p=getpos(a, i, j)) != NULL; j++)
			buf[p->colno] = p->v;
		for (j = 0; j < a->ncol; j++)
			printf("%10.6f ", buf[j]);
		printf("\n");
	}
	printf("\n");
}

void upper(Spm a, Elm *b) {
	Size i, j, k, basecol;
	struct spel *p, *q;
	Elm e;

	for (i = 0; i < a->nrow-1; i++)
		for (j = i+1; j < a->nrow; j++) {
			if ((p=find(a, j, i)) == NULL)
				continue;
			e = p->v / a->val[i].v;
			delete(a, j, i);
			for (k = 1; (p=getpos(a, i, k)) != NULL; k++) {
				basecol = p->colno;
				msg(stderr, "mtx: upper: doing %d, %d\n", i, basecol);
				if ((q=find(a, j, basecol)) != NULL)
					q->v -= p->v * e; 
				else
					add(a, j, basecol, -p->v * e);
			}
			msg(stderr, "mtx: upper: \n");
			b[j] -= b[i] * e;
		}
}

void lower(Spm a, Elm *b) {
	Size i, j, k, col;
	Elm e;
	struct spel *p, *q;

	for (i = a->nrow-1; i > 0; i--)
		for (j = i-1; j >= 0; j--) {
			if ((p=find(a, j, i)) == NULL)
				continue;
			e = p->v / a->val[i].v;
			delete(a, j, i);
			b[j] -= b[i] * e;
		}
}

void gauss_jordan(Spm a, Elm *b) {
	Size i;

	upper(a, b);
	lower(a, b);
	for (i=0; i < a->nrow; i++)
		b[i] /= a->val[i].v;
}

void printv(Elm *b, Size n) {
	while (n-- > 0)
		printf("%10.6f\n", *b++);
	printf("\n");
}

void sptom(Spm sm, Mtx m) {
	Size i, j;
	struct spel *p;

	for (i = 0; i < m->nrow; i++)
		for (j = 0; j < m->ncol; j++)
			if ((p=find(sm, i, j)) == NULL)
				m->val[i][j] = 0;
			else
				m->val[i][j] = p->v;
}

/* mtxread: read from element array val to matrix m */
void mtxread(Mtx m, Elm *val) {
	Size i, j;

	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < m->ncol; ++j)
			m->val[i][j] = *val++;
} 

/* mtxwrite: write matrix m to element array val */
void mtxwrite(Mtx m, Elm *val) {
	Rv *row;
	Rv p;
	
	for (row = m->val; row < m->val + m->nrow; ++row)
		for (p = *row; p < *row + m->ncol; ++p, ++val)
			*val = *p;
}

/* mtxcopy() - copy m2's val to m1 */
void mtxcopy(Mtx m1, Mtx m2) {
	Size i, j;

	for (i = 0; i < m1->nrow; ++i)
		for (j = 0; j < m1->ncol; ++j)
			m1->val[i][j] = m2->val[i][j]; /* a pointer */
}

Mtx mtxdup(Mtx a) {
	Mtx b = mtxalloc(a->nrow, a->ncol);

	if (b != NULL)
		mtxcopy(b, a);
	return b;
}

/* mtxmul -  return m1 * m2 */ 
Mtx mtxmul(Mtx m1, Mtx m2) {
	Size i, j, k;
	Mtx tmp;

	tmp = mtxalloc(m1->nrow, m2->ncol);
	for (i = 0; i < m1->nrow; ++i)
		for (j = 0; j < m2->ncol; ++j) {
			tmp->val[i][j] = 0;
			for (k = 0; k < m1->ncol; ++k)
				tmp->val[i][j]+=m1->val[i][k]*m2->val[k][j];
		}
	return tmp;
}

/* mtxdif: return m1 - m2 */
Mtx mtxdif(Mtx m1, Mtx m2) {
	Size i, j;
	Mtx tmp;

	tmp = mtxalloc(m1->nrow, m1->ncol);
	for (i = 0; i < m1->nrow; ++i)
		for (j = 0; j < m1->ncol; ++j)
			tmp->val[i][j] = m1->val[i][j] - m2->val[i][j];
	return tmp;
}

Mtx trans2(Mtx m) {
	Size i, j;
	Elm tmp;

	for (i = 0; i < m->nrow; ++i)
		for (j = i+1; j < m->ncol; ++j) {
			tmp = m->val[i][j];
			m->val[i][j] = m->val[j][i];
			m->val[j][i] = tmp;
		}
	return m;
}

Mtx trans(Mtx m) {
	Mtx a;

	a = mtxdup(m);
	if (a != NULL)
		trans2(a);
	return a;
}

/* return submatrix m(idx of row, :) */
Mtx rowsub(Mtx m, Size nrow, Size *idx) {
	Size i, j;
	Mtx tmp;

	tmp = mtxalloc(nrow, m->ncol);
	for (i = 0; i < nrow; ++i)
		for (j = 0; j < m->ncol; ++j)
			tmp->val[i][j] = m->val[idx[i]][j];
	return tmp;
}

void rowdup(Mtx sm, Mtx m, Size nrow, Size *in) {
	Size i, j;
	
	for (i = 0; i < nrow; ++i)
		for (j = 0; j < m->ncol; ++j)
			sm->val[i][j] = m->val[in[i]][j];
}

/* submatrix m(:, idx of col) */
Mtx colsub(Mtx m, Size ncol, Size *idx) {
	Size i, j;
	Mtx tmp;

	tmp = mtxalloc(m->nrow, ncol);
	for (i = 0; i < m->nrow; ++i)
		for (j = 0; j < ncol; ++j)
			tmp->val[i][j] = m->val[i][idx[j]];
	return tmp;
}

void coldup(Mtx sm, Mtx m, Size ncol, Size *in) {
	Size i, j;

	for (i = 0; i < m->nrow; i++)
		for (j = 0; j < ncol; j++)
			sm->val[i][j] = m->val[i][in[j]];
}

/* mtxprt() - print a matrix */
void mtxprint(char *title, Mtx a) {
	Size i, j;
	
	printf("%s =\n", title);
	for (i = 0; i < a->nrow; ++i) {
		for (j = 0; j < a->ncol; ++j)
			printf("%10.5f ", a->val[i][j]);
		putchar('\n');
	}
	putchar('\n');
}

void mtxnprt(int n, char *ttl, ...) {
	va_list ap;

	va_start(ap, ttl);
	while (n-- > 0)
		mtxprint(ttl, va_arg(ap, Mtx));
	va_end(ap);
} 

/* mtxalloc - allocate and return a nrow by ncol matrix.
	if there is not enough memory, then do cleaning,
	and return NULL */
Mtx mtxalloc(Size nrow, Size ncol) {
	Size i;
	Mtx tmp;

	if ((tmp = (Mtx) malloc(sizeof(struct mtx))) == NULL)
		return NULL;
	tmp->nrow = nrow;
	tmp->ncol = ncol;
	if ((tmp->val = (Rv *)malloc(sizeof(Elm *) * nrow))==NULL) {
		printf("mtxalloc: not enough memory\n");
		free(tmp);
		return NULL;
	}
	for (i = 0; i < nrow; ++i)
		if((tmp->val[i] = (Rv)malloc(sizeof(Elm) * ncol))==NULL){
			printf("mtxalloc: not enough memory\n");
			break;
		}
	if (i < nrow) {	/* not enough memory, do cleaning */
		while (--i >= 0)
			free(tmp->val[i]);
		free(tmp->val);
		free(tmp);
		return NULL;
	} 
	return tmp;
}

/* mtxfree() - free mtx m */ 
void mtxfree(Mtx m) {
	Size i;

	if (m == NULL)
		return;
	for (i = 0; i < m->nrow; ++i)
		free(m->val[i]);
	free(m->val);
	free(m);
}

/* exrow - exchange row p and q of matrix m */
void exrow(Mtx m, Size p, Size q) {
	swapt(Rv, m->val[p], m->val[q]);
}

void excol(Mtx m, Size i, Size j) {
	Size p;

	for (p = 0; p < m->nrow; ++p)
		swapt(Elm, m->val[p][i], m->val[p][j]);
}	

void readcol(Mtx m, Size colno, Elm *colv) {
	Size i;

	for (i = 0; i < m->nrow; ++i)
		m->val[i][colno] = colv[i];
}

int ludcmp(Mtx a, Size *indx, Elm *d) {
	Rv vv = (Rv) malloc(sizeof(Elm) * a->nrow);	
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

Mtx inverse(Mtx a) {
	Mtx b, inva;
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
	return trans2(inva);
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
