/* filename: 	pf.c 
   usage:	function definitions of Newton-Raphson power flow 
   changelog:
   2017-4-19	review
   */
#include <stdio.h>
#include <cs.h>
#include <klu.h>
#include "mtx.h"
#include "pwf.h"

/* makejac: create and update the Jacobian matrix */
cs *makejac(int dim) {
	static cs *jac_trip = NULL;
	int i, j;
	double v;

	if (dim == 0) /* Reuse the OLD Jacobian */
		return jac_trip;
	else {	/* delete and create the Jacobian matrix */
		if (jac_trip != NULL) 
			cs_spfree(jac_trip);
		if ((jac_trip = cs_spalloc(dim, dim, dim, 1, 1))==NULL)
			return NULL;
	}
	/* update the Jacobian */
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			if (getjac(&v, i, j))
				cs_entry(jac_trip, i, j, v); 
		}
	}
	return jac_trip;
}

int jacsolve(cs *matrix, double *arg) {
	cs *tmp;
	klu_common   comm;
	klu_symbolic *symb;
	klu_numeric  *nume;
	int n, *ap, *ai; 
	double *ax;
	int solstat;

	tmp = cs_compress(matrix);
	n = tmp->n;
	ap = tmp->p;
	ai = tmp->i;
	ax = tmp->x;
	klu_defaults(&comm);
	symb = klu_analyze(n, ap, ai, &comm);
	nume = klu_factor(ap, ai, ax, symb, &comm);	
	solstat = klu_solve(symb, nume, n, 1, arg, &comm);
	klu_free_symbolic(&symb, &comm);
	klu_free_numeric(&nume, &comm);
	return solstat;
}

void printjac(void) {
	printf("System Jacobian:\n");
	cs_print(makejac(0), 1);
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
	cs *jacob;

	flatstart(&errf);
	if (errf <= tol)
		return 0;
	dim = makeindex(&errf, &dfargs);
	for (i = 1; i <= lim; i++) {
		if ((jacob = makejac(dim)) == NULL)
			return -1;
		jacsolve(jacob, dfargs);
		updateindex(&errx);
		if (errx <= tol)
			break;
		updatenp();
		if (ischeck)
			checknode();
		dim = makeindex(&errf, &dfargs);	/* rebuild the index */
		if (errf <= tol)
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
