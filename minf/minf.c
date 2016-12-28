#include <stdio.h>
#include <stdlib.h>
#define ITLIM 100
#define ABS(a) ((a) >= 0 ? (a) : -(a)) 

typedef struct gen {
	double pg, min, max, a, b, c;
} gen;

int minf(int itlim, double eps, double step, 
	double pd, int ng, gen *sys, double *lamda, double *f);
void lmdtog(double lamda, gen *g);
double syspg(int n, gen *g);
double c(gen *g);
double ic(gen *g);
int inig(gen *g);

main () {
	gen *sys;
	int i, ng;
	double lamda, step, eps, f, pd, pdmax, pdmin;

	printf("input load range:\n");
	scanf("%lf-%lf\n", &pdmin, &pdmax);
	printf("input number of generators:\n");
	scanf("%d\n", &ng);
	sys = (gen *)malloc(sizeof(gen) * ng);
	printf("input generator data:\n");
	printf("pg\tpgmax\tpgmin\ta(pg^2)\tb(pg)\tc\n");
	for (i = 0; i < ng; ++i)
		inig(sys+i);
	printf("input step of lamda and eps:\n");
	scanf("%lf, %lf\n", &step, &eps);
	printf("load range: %7.2fMW - %7.2fMW\n", pdmin, pdmax);
	printf("num of gens: %d\n", ng);
	printf("step of lamda: %.3f$/MWh\n", step);
	printf("output power eps: %.3fMW\n", eps);

	printf("%7s %7s ", "lamda", "PD");
	for (i = 0; i < ng; ++i)
		printf("%5s%2d ", "PG", i+1);
	for (i = 0; i < ng; ++i)
		printf("%5s%2d ", "IC", i+1);
	printf("%7s\n", "fuel");
	printf("%7s %7s ", "$/MWh", "MW");
	for (i = 0; i < ng; ++i)
		printf("%7s ", "MW");
	for (i = 0; i < ng; ++i)
		printf("%7s ", "$/MWh");
	printf("%7s\n", "$/h");

	for (pd = pdmin; pd <= pdmax; pd += 10)
		if (minf(ITLIM, eps, step, pd, ng, sys, &lamda, &f) == 0) {
			printf("%7.2f %7.2f ", lamda, pd);
			for (i = 0; i < ng; ++i)
				printf("%7.2f ", sys[i].pg);
			for (i = 0; i < ng; ++i)
				printf("%7.2f ", ic(sys+i));
			printf("%7.2f", f);
			for (i = 0; i < ng && lamda == ic(sys+i); ++i)
				;
			if (i == ng)
				printf(" eco\n");
			else
				putchar('\n');
		} else
			printf("result is not convergent\n");
	return 0;	
}

int minf(int itlim, double eps, double step, 
	double pd, int ng, gen *sys, double *lamda, double *f) {

	int i, count;
	double spg;

	*lamda = ic(sys);
	for (i = 1; i < ng; ++i)
		if (*lamda > ic(sys+i))
			*lamda = ic(sys+i);

	for (count = 0; count < itlim; ++count) {
		spg = syspg(ng, sys);
		if (ABS(spg - pd) <= eps) {
			*f = 0;
			for (i = 0; i < ng; ++i)
				*f += c(sys+i);
			return 0;
		} else if (spg > pd)
			*lamda -= step;
		else
			*lamda += step;	
		for (i = 0; i < ng; ++i)
			lmdtog(*lamda, sys+i);
	}
	return -1;
}

void lmdtog(double lamda, gen *g) {
	double pg;

	pg = (lamda - g->b) / 2 / g->a;
	if (pg > g->max)
		pg = g->max;
	else if (pg < g->min)
		pg = g->min;
	g->pg = pg;
}

double ic(gen *g) {
	return 2 * g->a * g->pg + g->b;
}

double c(gen *g) {
	return g->a * g->pg * g->pg + g->b * g->pg + g->c;
}

double syspg(int n, gen *g) {
	int i;
	double pg;

	for (pg = i = 0; i < n; ++i)
		pg += g[i].pg;
	return pg;
}

int inig(gen *g) {
	scanf("%lf %lf %lf %lf %lf %lf\n",
	       &(g->pg), &(g->min), &(g->max), &(g->a), &(g->b), &(g->c));
	return 0;
}
