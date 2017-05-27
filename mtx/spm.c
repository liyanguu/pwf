#include <stdio.h> 
#include <stdlib.h>
#include "spm.h"
#include "msg.h"

/* spalloc - allocate a sparse matrix structure */
Spm spalloc(void) {
	return malloc(sizeof(struct spm_struct));
}

/* spelalloc - allocate a sparse element structure */
Spel spelalloc(void) {
	return malloc(sizeof(struct spel_struct));
}

/* spheadalloc - allocate a head of n elements */
Spel *spheadalloc(Size n) {
    	return malloc(sizeof(Spel) * n);
}

/* spinit - initialize an empty sparse matrix */
Spm spinit(Size nrow, Size ncol) {
	Spm sp;
	Size i;

	if ((sp = spalloc()) == NULL
	|| (sp->rowhead = spheadalloc(nrow)) == NULL
	|| (sp->colhead = spheadalloc(ncol)) == NULL)
	    	return NULL;
	sp->nrow = nrow;
	sp->ncol = ncol;
	for (i = 0; i < nrow; i++)
	    	sp->rowhead[i] = NULL;
	for (i = 0; i < ncol; i++)
	    	sp->colhead[i] = NULL;
	return sp;
}

Spm spinsert(Spel *head, Size i, Size j, Elm v) {
}
