#include <stdio.h>
#include "pwf.h"
#include "msg.h"

int main() {
	double ex, ef;

	print_test_msg = 1;

	readcdf(stdin);
	printlinef();
	printnode();
}
