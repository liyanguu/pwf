#include "simp.h"

int main() {
	LP lp;
	double cm[] = {
		3, 4, 1, 0, 0,
		5, -2, 0, 1, 0,
		2, 1, 0, 0, 1
	};
	double prof[] = { 16, 3, 0, 0, 0 };
	double solu[] = { 0, 0, 2, 1, 1 };
	int stat;

	lp = makelp(3, 5, cm, prof, solu);
	printlp(lp);
	while ((stat = simplex(lp)) && stat != UNBD)
		;
	if (!stat)
		printlp(lp);
	return 0;
}
