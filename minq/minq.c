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
double *q;
double s[100];
FILE *lf; 

void linearp(void);

int main() {
	Size maxnode, maxbranch, i;

	lf = fopen("../res/ieee14", "r");
	if (lf == NULL)
		return 1;
	readcdf(lf);
	pf(100, 1e-5, "n");
	getsize(&maxnode, &maxbranch);
	for (i = 0; i < maxnode; i++)
		normal_volt[i] = getnodeinfo(i, "volt");
	clear();
	closecdf(lf);
	
	lf = fopen("../res/ieee14.f4", "r");
	if (lf == NULL)
		return 2;
	readcdf(lf);
	linearp();
	return 0;
}

void linearp(void) {
	Size maxnode, maxbranch;
	lprec *lp;
	Mtx b, x;
	int n, i, j;
	int maxlim = 1;
	struct node *ti, *tj;
	
	getsize(&maxnode, &maxbranch);
	b = mtxalloc(maxnode, maxnode);
	for (n = 1; n <= maxlim; n++) {
		pf(100, 1e-5, "n");
		for (i = 0; (ti=getnode(i)) != NULL; i++) {
			ac_volt[i] = getnodeinfo(i, "volt");
			printf("ac volt: %f\n", ac_volt[i]);
			printf("ac volt: %f, %f\n", ti->volt.x, ti->volt.y);
			profarray[i+1] = 1;
		}
		for (i = 0; (ti=getnode(i)) != NULL; i++)
			for (j = 0; (tj=getnode(j)) != NULL; j++)
				b->val[i][j] = getsysinfo(ti, tj, 'l')/ac_volt[j];
		mtxprint("B", b);
		x = inverse(b);
		mtxprint("X", x);

		lp = make_lp(0, maxnode);
		set_obj_fn(lp, profarray);
		for (i = 0; i < maxnode; i++) 
			if ((vdiff[i]=normal_volt[i] - ac_volt[i]) > .0) {
				for (j = 1; j <= maxnode; j++)
					rowbuf[j]=x->val[i][j-1];
				add_constraint(lp, rowbuf, GE, vdiff[i]);
			}
		solve(lp);
		get_primal_solution(lp, s);
		totq = s[0];
		q = s + 1 + get_Nrows(lp);
		for (i = 0; i < maxnode; i++)
			setnodeinfo(i, q[i], "qgen");
		mtxfree(x);
		delete_lp(lp);
	}
	if (n <= maxlim) {
		printf("linear programming done\n");
		writecdf(stdout);
	}
}
