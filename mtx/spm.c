#include <stdio.h> 
#include <stdlib.h>
#include "spm.h"
#include "msg.h"

/* spalloc - allocate a sparse matrix structure */
Spm spalloc(void) {
	return (Spm) malloc(sizeof(struct spm_struct));
}

/* spelalloc - allocate a sparse element structure */
spel spelalloc(void) {
	return (spel) malloc(sizeof(struct spel_struct));
}

/* spcreat - create an empty sparse matrix */
Spm spcreat(Size nrow, Size ncol) {
	Spm sp;

	if ((sp = spalloc()) != NULL) {
		sp->elem = NULL;
		sp->nrow = nrow;
		sp->ncol = ncol;
		sp->nnz = 0;
		sp->cnz = (Size *) malloc(sizeof(Size) * ncol);
		sp->cnz_cnt = (Size *) malloc(sizeof(Size) * ncol);
		if (sp->cnz == NULL || sp->cnz_cnt == NULL)
			return NULL;
	}
	return sp;
}

/* spins - add a new element or update the elem's value */
int spins(Spm s, Size i, Size j, Elm v) {
	sprem(s, i, j);
	return spadd(s, i, j, v);
}

/* spadd - add an element (i, j, v) to s */
int spadd(Spm s, Size i, Size j, Elm v) {
	spel elem;

	if (!validate(s, i, j))
		return 0;
	if ((elem = spelalloc()) == NULL)
		return 0;
	elem->rowno = i;
	elem->colno = j;
	elem->v = v;
	elem->next = s->elem;
	s->elem = elem;
	++s->nnz;
	++s->cnz[j];
	return 1;
}

/* spfind - find and return the element ptr to [rowno, colno] of s */
spel  spfind(Spm s, Size rowno, Size colno) {
	spel  p;

	for (p=s->elem; p != NULL; p = p->next)
		if (p->rowno == rowno && p->colno == colno)
			break;
	return p;
}

/* spget - get the value of element [i, j] of s */
int spget(Spm s, Size i, Size j, Elm *v) {
	spel p;

	if ((p = spfind(s, i, j)) != NULL) {
		*v = p->v;
		return 1;
	}
	return 0;
}

/* sprem - remove an element [rowno, colno] of s */
void sprem(Spm s, Size rowno, Size colno) {
	spel p, pprev;

	for (pprev = NULL, p = s->elem; p != NULL; pprev = p, p = p->next)
		if (p->rowno == rowno && p->colno == colno) {
			if (pprev == NULL)
				s->elem = p->next;
			else
				pprev->next = p->next;
			free(p);
			s->nnz--;
			s->cnz[colno]--;
			return;
		}
	
}

/* spdelet - delete the sparse matrix s */
void spdelet(Spm s) {
	spel  p, q;

	if (s == NULL)
		return;
	for (p = s->elem; p != NULL; p = q) {
		printf("deleting %d, %d, %g\n", p->rowno, p->colno, p->v);
		q = p->next;
		free(p);
	}
	free(s);
}

void spprint(char *title, Spm s) {
	Size i, j;
	Elm v;

	if (s == NULL)
		return;
	printf("%s =\n", title);
	printf("non-zeros: %d\n", s->nnz);
	for (j = 0; j < s->ncol; j++)
		printf(" %10d", s->cnz[j]);
	printf("\n");
	for (i=0; i < s->nrow; i++) {
		for (j=0; j < s->ncol; j++) {
			printf(" %10.5f", spget(s, i, j, &v) ? v : 0);
		}
		printf("\n");
	}
	printf("\n");
}

void spcmp(Spm s, Size *nnz, Size *ap, Size *ai, Elm *ax) {
	spel p;
	Size j, tmp;

	*nnz = ap[s->ncol] = s->nnz;
	for (j = 0; j < s->ncol; j++) {
		s->cnz_cnt[j] = 0;
		ap[j] = 0;
		for (tmp = 0; tmp < j; tmp++) 
			ap[j] += s->cnz[tmp];
	}
	for (p = s->elem; p != NULL; p = p->next) {
		j = p->colno;
		tmp = ap[j] + s->cnz_cnt[j];	
		ai[tmp] = p->rowno;
		ax[tmp] = p->v;
		s->cnz_cnt[j]++;
	}
}
