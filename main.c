#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pwf.h"
#include "msg.h"

#define DEFLIM 150
#define DEFTOL 1e-3 
#define ERR_NO_INPUT_FILE_NAME "错误：没有指定输出文件名\n"
#define ERR_CANT_OPEN_FILE     "错误：无法打开输出文件 %s\n"
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
FILE *getfile(int argc, char **argv);

/* Gauss-Seidel / Newton power flow in IEEE common data format */
int main(int argc, char *argv[]) {
	char *progname = argv[0]; 
	char *p;
	int c;
	int stoploop = 0;

	ofp = stdout;
	lim = DEFLIM;
	tol = DEFTOL;
	while (--argc > 0) {
		if ((c = (*++argv)[0]) != '-' && c != '+')
			break;
		if (c == '+') {
			if (sscanf(&argv[0][1], "%lf", &tol) != 1)
				argc = 0;
			else if (tol < 0)
				tol = DEFTOL;
		} else if (sscanf(&argv[0][1], "%d", &lim) == 1) {
			if (lim < 1)
				lim = DEFLIM;
		} else {
			for (p = ++argv[0]; *p; ++p) {
				switch(*p) {
				case 't':
					print_test_msg = 1;
					if (!strcmp(++p, "pwf"))
						including("pwf");
					else if (!strcmp(p, "mtx"))
						including("mtx");
					stoploop = 1;
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
					ofp = getfile(--argc, ++argv);
					break;
				default:
					fprintf(stderr, ERR_WRONG_OPT, *p);
					argc = 0;
					stoploop = 1;
				}
				if (stoploop)
					break;
			}
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


FILE *getfile(int argc, char **argv) {
	FILE *ofp;

	if (argc <= 0) {
		fprintf(stderr, ERR_NO_INPUT_FILE_NAME);
		exit(3);
	} else if ((ofp = fopen(*argv, "w"))==NULL) {
		fprintf(stderr, ERR_CANT_OPEN_FILE, *argv);
		exit(3);
	}
	return ofp;
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
	} else if (flag.use_sparse) {
		meth = "s";
		fprintf(stderr, "稀疏牛顿法计算\n");
	} else if (flag.rect) {
		meth = "r";
		fprintf(stderr, "直角坐标牛顿法计算\n");
	} else {
		meth = "n";
		fprintf(stderr, "极坐标牛顿法计算\n");
	} 
	flatstart();
	inum = pf(lim, tol, meth, flag.adjusted);
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
