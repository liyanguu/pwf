#include "mtx.h"

typedef struct spxstru {
	lp_rec *lp;
	Mtx cm;
	Elm obj;
	Elm *cost;
	Elm *sol;
} Spx;

void makespx(Spx *px, Size nrow, Size ncol, Elm *array, Elm *cost) {
	int i;

	if ((px->lp=make_lp(nrow, ncol)) == NULL
	|| (px->cm=mtxalloc(nrow, ncol+1)) == NULL
	|| (px->cost=(Elm *) malloc(sizeof(Elm) * (ncol+1)))) == NULL
	|| (px->sol=(Elm *) malloc(sizeof(Elm) * (ncol+1)))) == NULL) 
		return;
	arraycpy(&px->cost[1], cost);
	for (i = 0; i < ncol; i++) {
		arraycpy(&px->cm->val[i][1], array);
		array += ncol;
	}
}
