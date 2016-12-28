#include <stdio.h>
#include <stdarg.h>
#include <string.h>

int print_test_msg = 0;

#define MAXNAME 100
static char *excludelist[MAXNAME];
static int pe = 0;

void exclude(char *name) {
	if (pe < MAXNAME)
		excludelist[pe++] = name;
	else
		pe = 0;
}

int inlist(char *name) {
	int i;

	for (i = 0; i < pe; i++)
		if (strstr(name, excludelist[i]) != NULL)
			return 1;
	return 0;
}

void msg(FILE *fp, char *fmt, ...) {
	va_list ap;

	if (!print_test_msg || inlist(fmt))
		return;
	va_start(ap, fmt);
	vfprintf(fp, fmt, ap);
	va_end(ap); 
}
