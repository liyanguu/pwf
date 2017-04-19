#include <stdio.h>
#include "lp_lib.h"

main(int argc, char *argv[]) {
	lprec *lp;
	FILE *fp;

	if ((fp = fopen("lptest.log", "a")) == NULL) {
		printf("can't open log\n");
		return 1;
	}
	fprintf(fp, "log info:\n");	
	lp = read_MPS(argv[1], NORMAL);
	if (lp == NULL) {
		fprintf(fp, "error: can't read model\n");
		return 1;
	}
	set_BFP(lp, "bfp_GLPK");
	switch(solve(lp)) {
	case 0:
		print_objective(lp);
		print_solution(lp, 4);
		break;
	case 2:
		printf("infeasible\n");
		return -2;
	default:
		break;
		return -1;
	}

	set_int(lp, 1, TRUE);
	write_LP(lp, fp);
	switch(solve(lp)) {
	case 0:
		print_objective(lp);
		print_solution(lp, 4);
		print_duals(lp);
		break;
	case 2:
		printf("infeasible\n");
		return -2;
	default:
		break;
		return -1;
	}

	fclose(fp);
	delete_lp(lp);
	return 0;
}
