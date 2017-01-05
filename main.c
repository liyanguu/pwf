#include <stdio.h>
#include <stdlib.h>
#include "pwf.h"
#include "msg.h"

#define DEFLIM 150
#define DEFERR 1e-3 

FILE *ifp, *ofp;
int lim;
double tol;
struct {
	unsigned int use_sparse : 1;
	unsigned int use_gauss : 1;
	unsigned int print_cdf : 1;
	unsigned int print_ybus : 1;
	unsigned int print_lineflow : 1;
	unsigned int print_node : 1;
} flag;

void loadflow(FILE *ifp, FILE *ofp);

/* gauss-seidel power flow in ieee common data format */
int main(int argc, char *argv[]) {
	char *progname = argv[0], *p;
	int c;

	ofp = stdout;
	lim = DEFLIM;
	tol = DEFERR;
	while (--argc > 0 && ((c = (*++argv)[0]) == '-' || c == '+'))
		if (c == '+') {
			if (sscanf(&argv[0][1], "%lf", &tol) != 1)
				argc = 0;
			else if (tol < 0)
				tol = DEFERR;
		} else if (sscanf(&argv[0][1], "%d", &lim) == 1) {
			if (lim < 1)
				lim = DEFLIM;
		} else {
			p = argv[0];
			while (c = *++p)
				switch(c) {
				case 't':
					print_test_msg = 1;
					if (strcmp(++p, "pwf") == 0)
						exclude("mtx");
					else if (strcmp(p, "mtx")==0)
						exclude("pwf");
					p = "a"; /* stop while */
					break;
				case 's':
					flag.use_sparse = 1;
					break;
				case 'g':
					flag.use_gauss = 1;
					break;
				case 'c':
					flag.print_cdf = 1;
					break;
				case 'y':
					flag.print_ybus = 1;
					break;
				case 'l':
					flag.print_lineflow = 1;
					break;
				case 'n':
					flag.print_node = 1;
					break;
				case 'o':
					if (--argc <= 0 ||
					(ofp = fopen(*++argv, "w"))==NULL) {
						fprintf(stderr, "错误：无法打开文件\n");
						exit(3);
					}
					break;
				default:
					argc = 0;
					fprintf(stderr, "错误：无效选项 %c\n", c);
					break;
				}
		}
	if (argc < 0) {
		fprintf(stderr, 
			"用法: %s -sgcnly -lim +tol -o outfile infile\n",
			progname); 
		exit(1);
	}
	if (argc == 0)
		loadflow(stdin, ofp);
	else 
		while (argc-- > 0) {
			if ((ifp = fopen(*argv, "r")) == NULL) {
				fprintf(stderr, 
					"error: can't open %s\n", *argv);
				exit(2);
			}
			argv++;
			loadflow(ifp, ofp);
			closecdf(ifp);
		}
	exit(0);
}

void loadflow(FILE *ifp, FILE *ofp) {
	double errf, errx;
	int inum;

	if (!readcdf(ifp))
		exit(2);
	fprintf(stderr, "迭代上限 %d, 容差 %10.8f\n", lim, tol);

	if (flag.use_gauss == 1)
		inum = pf(lim, tol, "g");
	else if (flag.use_sparse == 1)
		inum = pf(lim, tol, "s");
	else
		inum = pf(lim, tol, "n");

	if (inum >= 0 && inum <= lim) {
		fprintf(stderr, "迭代完成, 共%d次\n", inum);
	} else if (inum > lim)
		fprintf(stderr, "计算不收敛\n");
	else
		fprintf(stderr, "错误的计算方法\n");
	if (flag.print_cdf == 1)
		writecdf(ofp);
	if (flag.print_ybus == 1)
		printybus();
	if (flag.print_lineflow == 1)
		printlinef();
	if (flag.print_node == 1)
		printnode();
	clear();
}
