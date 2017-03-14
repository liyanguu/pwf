#include <stdio.h>
#include "pwf.h"

int main() {
	int it;

	readcdf(stdin);
	it = pf(100, 1e-4, "n", 0);
	reorder(NODETYPE);
	writecdf(stdout);
}
