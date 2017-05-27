/* main.c
   main function of pwf
   2017-4-25  review
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pwf.h"
#include "msg.h"

#define DEFLIM 100
#define DEFTOL 1e-3 
#define ERR_NO_INPUT_FILE "错误：没有指定输出文件\n"
#define ERR_CANT_OPEN_FILE     "错误：无法打开文件 %s\n"
#define ERR_WRONG_OPT          "错误：无效选项 %c\n"

FILE *ifp, *ofp;	/* load data input file and output file */
int lim;		/* iteration limit */
double tol;		/* tolerence */
struct {
	unsigned int use_sparse : 1;
	unsigned int use_gauss : 1;
	unsigned int print_cdf : 1;
	unsigned int print_ybus : 1;
	unsigned int print_lineflow : 1;
	unsigned int print_node : 1;
	unsigned int print_jac : 1;
	unsigned int adjusted : 1;
	unsigned int rect : 1;
} flag;

void loadflow(FILE *ifp, FILE *ofp);
int getopt(char *);

/* Gauss-Seidel / Newton power flow in IEEE common data format */
int main(int argc, char *argv[]) {
	char *progname = argv[0]; 
	char **infileptr;
	int c;
	int nfile;

	ofp = stdout;
	lim = DEFLIM;
	tol = DEFTOL;
	while (--argc > 0) {
		++argv;
		if ((c = (*argv)[0]) != '-' && c != '+')
			break;
		if (c == '+') {
			if (sscanf(&argv[0][1], "%lf", &tol) != 1)
				argc = 0;
			else if (tol < 0)
				tol = DEFTOL;
		} else if (c == '-') {
			if (sscanf(&argv[0][1], "%d", &lim) != 1) {
				if (!getopt(*argv))
					argc = 0;
			} else if (lim < 1)
				lim = DEFLIM;
		}
	}
	if (argc < 0) {
		fprintf(stderr, 
			"用法: %s [-agcnly] [-lim] [+tol] [-ooutfile] infile\n", progname); 
		exit(1);
	}

	nfile = argc;	/* processing the infiles */
	infileptr = argv;
	if (nfile == 0)
		loadflow(stdin, ofp);
	else 
		while (nfile-- > 0) {
			if ((ifp = fopen(*infileptr, "r")) == NULL) {
				fprintf(stderr, 
					ERR_CANT_OPEN_FILE, *argv);
				exit(1);
			}
			fprintf(stderr, "%s\n-----------\n", *infileptr);
			loadflow(ifp, ofp);
			closecdf(ifp);
			infileptr++;
		}
	exit(0);
}

void loadflow(FILE *ifp, FILE *ofp) {
	int inum;
	char *meth;

	if (!readcdf(ifp))
		exit(2);
	fprintf(stderr, "迭代上限 %d, 容差 %10.8f\n", lim, tol);

	if (flag.use_gauss) {
		meth = "g";
		fprintf(stderr, "高斯法计算\n");
	} else if (flag.rect) {
		meth = "r";
		fprintf(stderr, "直角坐标牛顿法计算\n");
	} else {
		meth = "n";
		fprintf(stderr, "极坐标牛顿法计算\n");
	} 

	inum = pf(meth, lim, tol, flag.adjusted);
	if (inum == -1) {
		fprintf(stderr, "错误的计算方法\n");
	} else if (inum > lim)
		fprintf(stderr, "计算不收敛\n");
	else
		fprintf(stderr, "迭代完成, 共%d次\n", inum);
		
	if (flag.print_cdf == 1)
		writecdf(ofp);
	if (flag.print_ybus == 1)
		printybus();
	if (flag.print_lineflow == 1)
		printlinef();
	if (flag.print_node == 1)
		printnode();
	if (flag.print_jac == 1)
		printjac();
	clear();
}


int getopt(char *p) {
	int c;

	while ((c = *++p) != '\0') {
		switch(c) {
		case 't':
			print_test_msg = 1;
			if (*++p)
				including(p);
			return 1;
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
		case 'j':
			flag.print_jac = 1;
			break;
		case 'a':
			flag.adjusted = 1;
			break;
		case 'r':
			flag.rect = 1;
			break;
		case 'o':
			if (*++p && (ofp = fopen(p, "w")) != NULL)
				return 1;
			if (*p == '\0')
				fprintf(stderr, ERR_NO_INPUT_FILE);
			else
				fprintf(stderr, ERR_CANT_OPEN_FILE, p);
			return 0; 	/* false */
		default:
			fprintf(stderr, ERR_WRONG_OPT, c);
			return 0;	/* false */
		}
	}
	return 1;	/* success */
}
