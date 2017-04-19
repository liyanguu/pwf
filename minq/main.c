#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pwf.h"
#include "mtx.h"
#include "msg.h"
#include "misc.h"

char *appendstr(char *base, char *append);
void exit_usg(char *progname);
char *getnextopt(int argc, char **argv);

int main(int argc, char **argv) {
	int c;
	int ncall, done;
	char *opt;
	char *casename = NULL;
	char *dataname = NULL;
	char *logname = NULL;
	char *progname = *argv;
	double tol = 0; 
	double alp = 0;
	
	while (--argc > 0) {
		if ((c=(*++argv)[0]) == '-' && (opt = getnextopt(argc, argv))!=NULL) {
			switch (c = *++argv[0]) {
			case 'l':
				logname = opt;
				break;
			case 'd':
				dataname = opt;
				break;
			case 't':
				tol = atof(opt);
				break;
			case 'a':
				alp = atof(opt);
				break;
			default:
				fprintf(stderr, "unknown option %c\n", c);
				argc = 0;
				break;
			}
			--argc;
			++argv;
		} else if (c != '-') {
			casename = *argv;
		} else {
			fprintf(stderr, "missing options\n");
			argc = 0;
		}
	}
	if (casename == NULL || argc < 0)
		exit_usg(progname);
	if (dataname == NULL)
		dataname = appendstr(casename, ".data");
	if (logname == NULL)
		logname = appendstr(casename, ".log");

	print_test_msg = 0;
	including("mtx");	/* print out mtx debug info */

	if (tol <= 0)
		tol = TOLR;
	minq_start(casename, logname, dataname);
	minq_makelp(minq_cnst());
	if (alp > 0)
		setalp(alp);
	done = 0;
	for (ncall = 1; ncall < 100 && !done; ncall++) {
		writelog("iteration: %d\n", ncall);
		if (!minq_lf())
			progend(-1, "load flow error\n");
		readvarin();
		done = minq_solve(tol);
		minq_log();
		writevarout();
	}
	if (done) {
		progend(0, "iteration done after %d times\n", ncall);	
	}
	progend(-1, "no solution\n");	
}


char *appendstr(char *base, char *append) {
	char *dst;

	dst = malloc(strlen(base) + strlen(append) + 1);
	if (dst != NULL) {
		strcpy(dst, base);
		strcat(dst, append);
	}
	return dst;
}

void exit_usg(char *progname) {
	fprintf(stderr, "usage: %s case [-l log] [-d data] [-t tol] [-a alp]\n",
		progname);
	exit(1);
}

char *getnextopt(int argc, char **argv) {
	return (--argc > 0) ? *++argv : NULL;
}
