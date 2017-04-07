#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pwf.h"
#include "mtx.h"
#include "msg.h"
#include "misc.h"

char *appendstr(char *base, char *append);

int main(int argc, char **argv) {
	int ncall, done;
	char *casename;
	char *dataname;
	char *logname;
	
	if (--argc > 0) {
		casename = *++argv;
		dataname = appendstr(casename, ".data");
		logname = appendstr(casename, ".log");
	}
	while (--argc > 0) {
		if ((*++argv)[0] == '-')
			switch (*++argv[0]) {
			case 'l':
				logname = ++argv[0];
				break;
			case 'd':
				dataname = ++argv[0];
				break;
			default:
				argc = 0;
				break;
			}
		else {
			progend(1, "usage: %s case [-llog] [-ddata]\n", *argv);
		}
	}
	if (argc < 0) {
		progend(1, "usage: %s case [-llog] [-ddata]\n", *argv);
	}

	print_test_msg = 1;
	including("mtx");	/* print out mtx debug info */

	minq_start(casename, logname, dataname);
	minq_makelp(minq_cnst());
	done = 0;
	for (ncall = 1; ncall < 100 && !done; ncall++) {
		writelog("iteration: %d\n", ncall);
		if (!minq_lf())
			progend(-1, "load flow error\n");
		readvarin();
		done = minq_solve(TOLR);
		minq_log();
		writevarout();
	}
	if (done) {
		progend(0, "iteration done\n");	
	}
	progend(-1, "no solution\n");	
}


char *appendstr(char *base, char *append) {
	char *dst;

	dst = malloc(strlen(base) + strlen(append) + 1);
	if (dst) {
		strcpy(dst, base);
		strcat(dst, append);
	}
	return dst;
}
