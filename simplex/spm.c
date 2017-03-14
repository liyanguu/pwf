#include <stdio.h> 
#include <stdlib.h>
#include "spm.h"

/* spalloc - allocate a sparse matrix structure */
Spm spalloc(void) {
	return (Spm) malloc(sizeof(struct spm));
}

/* spelalloc - allocate n element */
struct spel *spelalloc(Size n) {
	return (struct spel *) malloc(sizeof(struct spel) * n);
}

/* spcreat - create a nrow by ncol matrix */
Spm spcreat(Size nrow, Size ncol) {
	Spm m;
	Size i;

	m = spalloc();
	m->val = spelalloc(nrow);
	m->nrow = nrow;
	m->ncol = ncol;
	for (i=0; i < nrow; i++) {
		m->val[i].v = .0;
		m->val[i].colno = 0;
		m->val[i].nextcol = NULL;
	}
	return m;
}

/* sprow - return the pointer to the first none diagnal element */
struct spel *sprow(Spm a, Size rowno) {
	if (rowno >= 0 && rowno < a->nrow)
		return a->val[rowno].nextcol;
	else {
		fprintf(stderr, "sprow: rowno out of range %d\n", rowno);
		return NULL;
	}
}

/* sppos - get the element at position pos behind head */
struct spel *sppos(struct spel *head, Size pos) { 
	struct spel *p;

	for (p = head; pos-- > 0 && p != NULL; p=p->nextcol)
		;
	return p;
}

/* spfind - find the none-diagnal element (rowno, colno) */
struct spel *spfind(Spm a, Size rowno, Size colno) {
	struct spel *p;

	for (p=sprow(a, rowno); p != NULL; p = p->next)
		if (p->colno == colno)
			break;
	return p;
}

/* spinsert - insert or update element v at (rowno, colno) of a */
int spinsert(Spm a, Size rowno, Size colno, Elm v) {
	struct spel *p;

	if (rowno < 0 || rowno >= a->nrow ||
		colno < 0 || colno >= a->ncol || 
		v == 0)
		return -1;
	if (rowno == colno)		/* diagnal element */
		a->val[rowno].v = v;
	else if ((p = spfind(a, rowno, colno)) != NULL)	/* none diagnal */
		p->v = v;		/* if exists, replace it */
	else {
		p = spelalloc(1);
		p->colno = colno;
		p->v = v;
		p->nextcol = a->val[rowno].nextcol;
		a->val[rowno].nextcol = p;
		a->val[rowno].colno++;
	}
	return 0;
}

/* spremove - remove element (rowno, colno) */
void spremove(Spm a, Size rowno, Size colno) {
	struct spel *p, *q;

	if (rowno < 0 || rowno >= a->nrow ||
	  colno < 0 || colno >= a->ncol ||
	  colno == rowno)
		return;
	for (p = &a->val[rowno]; p->nextcol != NULL; p = p->nextcol)
		if (p->nextcol->colno == colno) {
			q = p->nextcol;
			p->nextcol = q->nextcol;
			fprintf(stderr, "spremove: deleting %d, %d\n", rowno,q->colno);
			free(q);
			a->val[rowno].colno--;
			return;
		}
}

/* spdeletrow - delete a row */
void spdeletrow(struct spel *p) {
	struct spel *q;

	for ( ; p != NULL; p = q) {
		q = p->nextcol;
		free(p);
	}
}

/* spdelet - delete a maritx */
void spdelet(Spm a) {
	Size i; 

	if (a == NULL)
		return;
	for (i = 0; i < a->nrow; i++)
		spdeletrow(sprow(a, i));	/* free none diagnal */ 
	free(a->val);	/* free diagnal */
	free(a);	/* free the structure */
}

void spprint(char *title, Spm a) {
	struct spel *p;
	Size i, j;

	printf("%s =\n", title);
	for (i=0; i < a->nrow; i++) {
		for (j=0; j < a->ncol; j++) {
			if (i == j)
				printf("%10.6f ", a->val[i].v);
			else if ((p = spfind(a, i, j)) != NULL)
				printf("%10.6f ", p->v);
			else
				printf("%10.6f ", 0.);
		}
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
	Size i, j;
	Elm e;
	struct spel *p;

	for (i = a->nrow-1; i > 0; i--) {
		for (j = i-1; j >= 0; j--) {
			if ((p=find(a, j, i)) == NULL)
				continue;
			e = p->v / a->val[i].v;
			delete(a, j, i);
			b[j] -= b[i] * e;
		}
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
		printf("%10.6f ", *b++);
	printf("\n");
}

void sptom(Spm sm, Mtx m) {
	Size i, j;
	struct spel *p;

	for (i = 0; i < m->nrow; i++)
		for (j = 0; j < m->ncol; j++)
			if ((p=spfind(sm, i, j)) == NULL)
				m->val[i][j] = 0;
			else
				m->val[i][j] = p->v;
}
