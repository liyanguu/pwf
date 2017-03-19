#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <fenv.h>

int print_test_msg = 0;

#define MAXLIST 100
static char *excludelist[MAXLIST+1];
static int pe = 0;
static char *includelist[MAXLIST+1];
static int pi = 0;

void excluding(char *name) {
	if (pe < MAXLIST)
		excludelist[pe++] = strdup(name);
}

void including(char *name) {
	if (pi < MAXLIST)
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
	#pragma STDC FENV_ACCESS on
	int res = fetestexcept(FE_ALL_EXCEPT);

	if (res & FE_INEXACT)
		msg(stderr,"warning: floating point: loss precision\n");
	if (res & FE_UNDERFLOW)
		msg(stderr,"warning: floating point: underflow\n");
	if (res & FE_OVERFLOW)
		msg(stderr,"warning: floating point: overflow\n");
	if (res & FE_DIVBYZERO)
		msg(stderr,"warning: floating point: div by zero\n");
	if (res & FE_INVALID)
		msg(stderr,"warning: floating point: invalid operation\n");
	feclearexcept(FE_ALL_EXCEPT);
}
