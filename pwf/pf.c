/* filename: 	pf.c 
   usage:	function definitions of Newton-Raphson power flow 
   changelog:
   2017-4-19	review
   2017-5-13    
   	6-1	colspm 版本
   */
#include <stdio.h>
#include "pwf.h"
#include "wrklu.h"
#include "colspm.h"

static ColSpm *jac_matrix = NULL;

/* makejac: create and update the Jacobian matrix */
ColSpm *makejac(int dim) {
	int i, j;
	double v;

	if (dim <= 0) /* Reuse the OLD Jacobian */
		return jac_matrix;
	else if (jac_matrix != NULL)
		colspm_clr(jac_matrix); /* refresh Jacobian matrix */
	else 
	    	jac_matrix = colspm_init(dim); /* build the Jacobian */
	/* update the Jacobian */
	for (j = 0; j < dim; j++) {
	    	colspm_addcol(jac_matrix, j);
		for (i = 0; i < dim; i++)
			if (getjac(&v, i, j))
			    	colspm_add(jac_matrix, i, v);
		colspm_endcol(jac_matrix);
	}
	colspm_adjust(jac_matrix);
	return jac_matrix;
}

void deletjac(void) {
    	jac_matrix = colspm_del(jac_matrix);
}

int jacsolve(ColSpm *jmatrix, double *arg) {
	int n, *ap, *ai; 
	double *ax;
	int solstat;

	n = jmatrix->dim;
	ap = jmatrix->pos;
	ai = jmatrix->rows;
	ax = jmatrix->elms;
	wr_klu_defaults();
	wr_klu_analyze(n, ap, ai);
	wr_klu_factor(ap, ai, ax);
	solstat = wr_klu_solve(n, 1, arg);
	return solstat;
}

void printjac(void) {
    	ColSpm *jmatrix = makejac(-1);

	if (jmatrix == NULL)
	    	printf("NULL Jacobian\n");
	else {
		printf("System Jacobian:\n");
		colspm_print(jmatrix);
	}
}

/* static int *indx;  swap index of ludcmp & lubksb */

/* nrpf: Newton-Raphson power flow
	lim: interation limit
	tol: tolerance
	ischeck:
		1  regulated load flow
		0  unregulated load flow
   return values:
	-1	error
	<= lim	iteration convergent
	> lim	iteration nonconvergent
*/
int nrpf(int lim, double tol, int ischeck) {
	/*
	dfargs: mismatch vetcor [DP, DQ] 
	errf: abs max of args
	errx: abs max of dA, dV or de , df
	*/
	int i, dim;
	double errf, errx;
	double *dfargs;

	flatstart(&errf);
	if (errf <= tol)
		return 0;
	for (i = 1; i <= lim; i++) {
		dim = makeindex(&errf, &dfargs);
		if (errf <= tol)
			break;
		if (makejac(dim) == NULL)
			return -1;
		jacsolve(makejac(-1), dfargs);
		updateindex(&errx);
		if (ischeck && checknode() > 0) /* rebuild the index */
		    	continue;
		if (errx <= tol)
			break;
	}
	return i;
}

int gspf(int lim, double tol, int ischeck) {
	int i;
	double errf;

	for (i = 1; i <= lim; i++) {
		if (!gs(&errf))
			return -1;
		if (ischeck && checknode() > 0)
			continue;
		if (errf <= tol)
			break;
	}
	return i;
}

/* pf: main function for power flow calculations
	method - 
		'n' - polar Newton-Raphson
		'r' - rectangular Newton-Raphson (not supported now)
		'g' - Gauss-Seidal
   return values:
   	same as nrpf.
*/

int pftype;

int pf(char *method, int lim, double tol, int ischeck) {
	switch(*method) {
	case 'n':
		pftype = NR_POL;
		return nrpf(lim, tol, ischeck);	
	case 'r':
		pftype = NR_REC;
		return nrpf(lim, tol, ischeck);	
	case 'g':
		return gspf(lim, tol, ischeck);
	default:
		return -1;
	}
}
