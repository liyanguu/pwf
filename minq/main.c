#include <stdio.h>
#include "pwf.h"
#include "mtx.h"
#include "msg.h"
#include "misc.h"

int main(int argc, char **argv) {
	int ncall, done;
	char *casename, *logname, *dataname;
	
	if (argc != 4) {
		progend(-1, "usage: %s case log data\n", *argv);
	}
	casename = argv[1];
	logname = argv[2];
	dataname = argv[3];

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
