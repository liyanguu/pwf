/* spm.h - SParse Matrix definition */

#ifndef header_spm
#define header_spm

#include "mtx.h"

typedef struct spel_struct {		/* a sparse matrix element */
	Elm v;		/* value */
	Size colno;	/* coloumn no */
	Size rowno;	/* row no */
	struct spel_struct *next;	/* ptr to the next element */
} *spel; 

typedef struct spm_struct {		/* a sparse matirx */
	spel elem;		/* ptr to elements */ 
	Size nrow;	/* # rows */
	Size ncol;	/* # coloumns */
	Size nnz;	/* # non-zero elements */
	Size *cnz;	/* # non-zero elements in a column */
	Size *cnz_cnt;	/* a counter for cnz */
	Size *ap;	/* compressed-column data */
	Size *ai;
	Elm *ax;
} *Spm;

#define validate(s, i, j) ((i) >= 0 && (i) < s->nrow && (j) >=0 && (j) < s->ncol)

void spprint(char *title, Spm a);
void spdelet(Spm a);
void sprem(Spm a, Size rowno, Size colno);
int spadd(Spm a, Size rowno, Size colno, Elm v);
int spget(Spm a, Size rowno, Size colno, Elm *v);
spel spfind(Spm a, Size rowno, Size colno);
Spm spcreat(Size nrow, Size ncol);
Spm spalloc(void);
spel spelalloc(void);
void spcmp(Spm s, Size *nz, Size *ap, Size *ai, Elm *ax);

#endif
