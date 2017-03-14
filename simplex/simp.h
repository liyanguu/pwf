/* simplex algo version 2 */
/* simp.h */

#ifndef H_SIMP
#define H_SIMP

#include "mtx.h"

enum { OPTM, UNDT, UNBD }; 
typedef struct lpstru {
	Mtx cm; 	/* constraint matrix */
	Mtx bm;         /* base matrix */
	Mtx bm_in;	/* base inverse matrix */
	Mtx nbm; 	/* none-base matrix */
	Mtx pv; 	/* profit vector */
	Mtx sv; 	/* solution vector */
	Mtx bsv;  	/* base solution vector*/
	Mtx bpv;        /* base profit vector */
	Mtx npv;	/* none-base profit vector */
	Size *idx;
	Size nrow, ncol;
	Mtx obj;
	int flag;
	/* calculation intermediate values: */
	Mtx dlt_tr;	/* delta value's transpose */	
	Mtx yv;		/* vector y */
	Mtx tnbm, tbm, tcm;
} *LP;

#define isopt(lp) ((lp)->flag == OPT)
#define isunbound(lp) ((lp)->flag == UNBD)
#define object(lp) ((lp)->obj->val[0][0])

LP makelp(Size nrow, Size ncol, Elm *constraint, Elm *profitco, Elm *solution);
void printlp(LP lp);
int simplex(LP lp);

#endif
