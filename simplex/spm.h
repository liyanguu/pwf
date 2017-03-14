/* spm.h - sparse matrix definition */

#ifndef H_SPARSE
#define H_SPARSE

#include "mtx.h"

struct spel {		/* a sparse matrix element */
	Elm v;		/* value */
	Size colno;	/* coloumn no */
	struct spel *nextcol;
};

typedef struct {		/* a sparse matirx */
	struct spel *val;	/* diagnal elements */ 
	Size nrow, ncol;	/* dimension */
} *Spm;

void gauss_jordan(Spm, Elm *);
void spprint(char *title, Spm a);
void spdelet(Spm a);
void spdeletrow(struct spel *p);
void spremove(Spm a, Size rowno, Size colno);
int spinsert(Spm a, Size rowno, Size colno, Elm v);
struct spel *spfind(Spm a, Size rowno, Size colno);
struct spel *sppos(struct spel *head, Size pos);
struct spel *sprow(Spm a, Size rowno);
Spm spcreat(Size nrow, Size ncol);
Spm spalloc(void);

#endif
