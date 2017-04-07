/* pf.c - implementation of Newton-Raphson power flow */

#include <stdio.h>
#include <cs.h>
#include <klu.h>
#include "mtx.h"
#include "pwf.h"

/* makejac - create and / or update the Jacobian matrix */
cs *makejac(int dim) {
	int i, j;
	double v;
	static cs *jac_trip;

	if (dim == 0) /* Reuse the OLD Jacobian */
		return jac_trip;
	else {	/* delete the old Jacobian matrix */
		if (jac_trip != NULL) 
			cs_spfree(jac_trip);
		if ((jac_trip = cs_spalloc(dim, dim, dim, 1, 1))==NULL)
			return NULL;
	}

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
	for (i = 1; i <= lim; i++) {
		if ((jacob = makejac(dim)) == NULL)
			return -1;
		jacsolve(jacob, dfargs);
		norm(dfargs, myabs(dim), &errx);
		if (errx <= tol)
			return i;
		updateindex();
		updatenp();
		if (ischeck)
			checknode();
		dim = makeindex(&errf, &dfargs);	/* rebuild the index */
		if (errf <= tol)
			return i;
	}
	return i;
}

int gspf(int lim, double tol, int ischeck) {
	int i;
	double errf;

	for (i = 1; i <= lim; i++) {
		if (!gs(&errf))
			return -1;
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
