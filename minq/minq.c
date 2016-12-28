#include <stdio.h>
#include "pwf.h"
#include "simplex/mtx.h"
#include "lpsolve55/lp_lib.h" 

double normal_volt[100];
double ac_volt[100];
double profarray[100];
double vdiff[100];
double rowbuf[100];
double totq;
double *q, *s;
FILE *lf; 

void linearp(void);

int main() {
	Size maxnode, maxbranch, i;

	lf = fopen("../resourse/ieee14", "r");
	if (lf == NULL)
		return 1;
	readcdf(lf);
	pf(100, 1e-5, "n");
	getsize(&maxnode, &maxbranch);
	for (i = 0; i < maxnode; i++)
		normal_volt[i] = getnodeinfo(i, "volt");
	clear();
	closecdf(lf);
	
	lf = fopen("../resourse/ieee14.f4", "r");
	if (lf == NULL)
		return 2;
	readcdf(lf);
	pf(100, 1e-5, "n");
	linearp();
	return 0;
}

void linearp(void) {
	Size maxnode, maxbranch;
	lprec *lp;
	Mtx b, x;
	int n, i, j;
	int maxlim = 1;
	
	getsize(&maxnode, &maxbranch);
	b = mtxalloc(maxnode, maxnode);
	for (n = 1; n <= maxlim; n++) {
		pf(100, 1e-5, "n");
		for (i = 0; i < maxnode; i++) {
			ac_volt[i] = getnodeinfo(i, "volt");
			profarray[i+1] = 1.;
		}
		for (i = 0; i < maxnode; i++)
			for (j = 0; j < maxnode; j++)
				b->val[i][j] = getsen(i, j, "l")/ac_volt[j];
		mtxprint("B", b);
		x = inverse(b);
		mtxprint("X", x);

		lp = make_lp(0, maxnode);
		set_obj_fn(lp, profarray);
		for (i = 0; i < maxnode; i++)
			if ((vdiff[i]=normal_volt[i] - ac_volt[i]) > .0) {
				for (j = 1; j <= maxnode; j++)
					rowbuf[j]=x->val[i][j-1];
				add_contraint(lp, rowbuf, GE, vdiff[i]);
			}
		solve(lp);
		get_primal_solution(lp, s);
		totq = *s;
		q = s + 1 + get_Nrows(lp);
		for (i = 0; i < maxnode; i++)
			setnodeinfo(i, q[i], "qgen");
		mtxfree(x);
		delet_lp(lp);
	}
	if (n <= maxlim) {
		printf("linear programming done\n"
		writecdf(stdout);
	}
}
