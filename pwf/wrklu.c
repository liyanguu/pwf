/* wrklu.c - KLU 函数的皮函数 
 * 2017-5-27 */

#include "wrklu.h"

static klu_common common;
static klu_symbolic *symbolic;
static klu_numeric *numeric;
static int dimension = 0;
static int set;

void wr_klu_defaults(void) {
    	if (set == 0) {
    		klu_defaults(&common);
		symbolic = NULL;
		numeric = NULL;
	    	set = 1;
	}
}

void wr_klu_analyze(int n, int *ap, int *ai) {
	if (dimension == n)
	    	return;
	if (dimension != n && dimension != 0)
		klu_free_symbolic(&symbolic, &common);
	symbolic = klu_analyze(n, ap, ai, &common);
	dimension = n;
}

void wr_klu_factor(int *ap, int *ai, double *ax) {
    	if (dimension > 0)
	    	klu_free_numeric(&numeric, &common);
	numeric = klu_factor(ap, ai, ax, symbolic, &common);
}

int wr_klu_solve(int nldim, int nrhs, double *sol) {
	return klu_solve(symbolic, numeric, nldim, nrhs, sol, &common);
}
