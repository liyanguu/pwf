/* pf.c - implementation of Newton-Raphson power flow */

#include <stdio.h>
#include "cs.h"
#include "klu.h"
#include "pwf.h"

cs *makejac(int dim) {
	int i, j;
	static cs *jac_trip;
	static int jac_dim;


	if (dim <= 0) /* use the OLD Jacobian */
		dim = jac_dim;
	else {/* delete the old Jacobian matrix */
		if (jac_trip != NULL) 
			cs_spfree(jac_trip);
		jac_trip = cs_spalloc(dim, dim, dim, 1, 1);
		jac_dim = dim;
	}

	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			cs_entry(jac_trip, i, j, getjac(i, j));
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
	dfargs -  mismatch vetcor [DP, DQ] 
	errf - max of args
	errx - max of dA, dV or de , df
   return values:
	>= 0 - success
*/
int nrpf(int lim, double tol, int ischeck) {
	int i, dim;
	double errf, errx;
	double *dfargs;
	cs *jacob;

	flatstart(&errf);
	if (errf <= tol)
		return 0;
	dim = makeindex(&errf, &dfargs);
	if (errf <= tol)
		return 0;
	for (i = 1; i <= lim; i++) {
		jacob = makejac(dim);
		jacsolve(jacob, dfargs);
		norm(dfargs, abs(dim), &errx); 
		if (errx <= tol)
			return i;
		updateindex();
		updatenp();
		if (ischeck && checknode() > 0)
			dim = makeindex(&errf, &dfargs);
		else
			dim = -dim;
	}
	return i;
}

int gspf(int lim, double tol, int ischeck) {
	int i;
	double errf;

	for (i = 1; i <= lim; i++) {
		gs(&errf);
		if (errf <= tol)
			break;
	}
	return i;
}

/* pf: power flow calculations
	lim - interation limit
	tol - tolerance
	method - 
		'n' - polar Newton-Raphson
		'r' - rectangular Newton-Raphson
		'g' - Gauss-Seidal
	ischeck -
		1 - regulated load flow
		0 - unregulated load flow
   return values:
	<= lim - convergence
	>  lim - none convergence
	-1 - errors
*/

int pftype;

int pf(int lim, double tol, char *method, int ischeck) {
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
