#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <fenv.h>

int print_test_msg = 0;

#define MAXNAME 100
static char *excludelist[MAXNAME+1];
static int pe = 0;
static char *includelist[MAXNAME+1];
static int pi = 0;

void excluding(char *name) {
	if (pe < MAXNAME)
		excludelist[pe++] = strdup(name);
}

void including(char *name) {
	if (pi < MAXNAME)
		includelist[pi++] = strdup(name);
}

int inlist(char *name, char **list) {
	char *p;

	while ((p = *list++) != NULL)
		if (strstr(name, p) != NULL)
			return 1;
	return 0;
}

void msg(FILE *fp, char *fmt, ...) {
	va_list ap;

	if (print_test_msg 
	&& (!inlist(fmt, excludelist) && inlist(fmt, includelist))) {
		va_start(ap, fmt);
		vfprintf(fp, fmt, ap);
		va_end(ap); 
	}
}

void print_except(void) {
	int res = fetestexcept(FE_ALL_EXCEPT);

	if (res & FE_INEXACT)
		printf("warning: floating point: loss precision\n");
	if (res & FE_UNDERFLOW)
		printf("warning: floating point: underflow\n");
	if (res & FE_OVERFLOW)
		printf("warning: floating point: overflow\n");
	if (res & FE_DIVBYZERO)
		printf("warning: floating point: div by zero\n");
	if (res & FE_INVALID)
		printf("warning: floating point: invalid operation\n");
	feclearexcept(FE_ALL_EXCEPT);
}
