/* spm.h - SParse Matrix definition */

#ifndef H_spm_H
#define H_spm_H

#include "mtx.h"

typedef struct spel_struct {		/* a sparse matrix element */
	Elm v;		/* value */
	Size colno;	/* coloumn no */
	struct spel_struct *nextcol;	/* ptr to the next element in a row */
} *spel; 

typedef struct spm_struct {		/* a sparse matirx */
	spel val;		/* diagnal elements */ 
	Size nrow, ncol;	/* dimension */
} *Spm;

#define validate(a, i, j) ((i) >= 0 && (i) < (a)->nrow && (j) >= 0 && (j) < (a)->ncol)

void gauss_jordan(Spm, Elm *);
void spprint(char *title, Spm a);
void spdelet(Spm a);
void spdeletrow(spel p);
void spremove(Spm a, Size rowno, Size colno);
int spinsert(Spm a, Size rowno, Size colno, Elm v);
spel spfind(Spm a, Size rowno, Size colno);
spel sppos(Spm a, Size rowno, Size pos);
spel sprow(Spm a, Size rowno);
Spm spcreat(Size nrow, Size ncol);
Spm spalloc(void);

#endif
