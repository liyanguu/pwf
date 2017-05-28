/* wrapper.c - KLU 函数的皮函数 
 * 2017-5-27 */

#include "wrapper.h"

static klu_common common;
static int set;
static klu_symbolic *symbolic;
static klu_numeric *numeric;

void wr_klu_defaults(void) {
    	if (set == 0) {
    		klu_defaults(&common);
	    	set = 1;
	}
}

void wr_klu_analyze(int n, int *ap, int *ai) {
    	static int dimension = 0;

	if (n <= 0) {
	    	symbolic = NULL;
	    	return;
	}
	if (dimension == n)
	    	return;
	if (dimension != n && dimension != 0)
		klu_free_symbolic(&symbolic, &common);
	symbolic = klu_analyze(n, ap, ai, &common);
	dimension = n;
}

void wr_klu_factor(int *ap, int *ai, double *ax) {
    	if (numeric != NULL)
	    	klu_free_numeric(&numeric, &common);
	numeric = klu_factor(ap, ai, ax, symbolic, &common);
}

int wr_klu_solve(int nldim, int nrhs, double *sol) {
	return klu_solve(symbolic, numeric, nldim, nrhs, sol, &common);
}
