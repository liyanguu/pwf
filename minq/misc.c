/* misc.c - minq auxiliary functions */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdarg.h>
#include <lp_lib.h>
#define LPSOLVEAPIFROMLIBDEF
#include <lp_explicit.h>
#include "mtx.h"
#include "pwf.h"
#include "misc.h"

FILE *logfp;
FILE *setfp;
FILE *lffp;
double vmin[MAXVAR];
double vmax[MAXVAR];
double qmax[MAXVAR];	/* qgen limits */
double theta[MAXVAR];	/* power factor limits */
double alp;		/* iteration factor alpha, < 1 */

REAL vcormin[MAXVAR];
REAL vcormax[MAXVAR];
REAL cm[MAXVAR];
REAL cost[MAXVAR];		/* LP 输入：费用系数 */
REAL vc[MAXVAR];        /* LP 结果：预计的电压补偿量 */
REAL qgen[MAXVAR];         /* LP 结果：节点无功出力向量 */
REAL totq;				/* LP 结果：总无功出力 */
REAL vac[MAXVAR];		/* 输入数据，潮流计算的电压值 */
REAL vdiff[MAXVAR];		/* 潮流计算电压与最小电压之偏差 */
REAL maxdiff;			/* 最大电压偏差 */
static int busno[MAXVAR];	/* 实际控制节点编号 (psasp 的内部编号)*/
static int dim;
lprec *lp;
hlpsolve lpsolve;	/* ptr to the lpsolve shared library */
#define LP_LIB "lpsolve55.so"	/* name of the lpsolve shared lib */

void readvarin(void) {
	int i;
	struct node *t;

	for (i = 0; i < dim; i++) {
		t = getnode(busno[i]);
		vac[i] = getnodeinfo(t, VOLT);
	}
}

void writevarout(void) {
	int i;

	for (i = 0; i < dim; i++) {
		setnodeinfo(getnode(busno[i]), QGEN, qgen[i]);
	}
}

void minq_start(char *casename, char *logname, char *dataname) {
		int i;

		if ((lffp = fopen(casename, "r")) == NULL)
			progend(1, "can't read load flow data %s\n", casename);
		if ((logfp = fopen(logname, "w")) == NULL)  /* 打开 log 文件 */
			progend(1, "can't open log file %s\n", logname);
		if ((setfp = fopen(dataname, "r")) == NULL) /* 打开设定数据文件*/
			progend(1, "can't open data file %s\n", dataname);
		if (!readdata(setfp))
			progend(1, "can't read data file %s\n", dataname);
		fclose(setfp);
		writelog("#Buses\t%d\n", dim);
		writelog("Busno\t");
		for (i = 0; i < dim; ++i)
			writelog("%10d ", busno[i]+1);
		writelog("\n");
		writelog("Vmin\t");
		writearray(logfp, vmin, dim);
		writelog("Vmax\t");
		writearray(logfp, vmax, dim);
		writelog("Qmax\t");
		writearray(logfp, qmax, dim);
		readcdf(lffp);
		minq_lf();
}

int minq_lf(void) {
	int it;

	if ((it = pf(LIMT, TOLR, "n", 0)) < 0) {
		writelog("load flow error\n");
		return 0;
	} else if (it > LIMT) {
		writelog("load flow none convergent\n");
		return 0;
	}
	return 1;
}

/* sysmatrix: 从输入向量读取系统的节点导纳矩阵。
  返回系统的节点B阵，并输出其与G阵 minq.matrix
*/
Mtx sysmatrix(void) {
	int i, j;
	int nbus;
	double v;
	Mtx mtx_b;
	FILE *fp;

	if ((fp = fopen("minq.matrix", "w")) == NULL) {
		writelog("can't open minq.matrix\n");
		return NULL;
	}
	if ((nbus = getnnode()) == 0) {
		writelog("no load flow result\n");
		return NULL;
	}
	fprintf(fp, "G part:\n");
	for (i = 0; i < nbus; ++i) {
		for (j = 0; j < nbus; ++j)
			if (getsysinfo(&v, getnode(i), getnode(j), G))
				fprintf(fp, "%10.5f ", v);
			else
				fprintf(fp, "%10s ", " ");
		fputc('\n', fp);
	}

	mtx_b = mtxalloc(nbus, nbus);
	fprintf(fp, "\nB part:\n");
	for (i = 0; i < nbus; ++i) {
		for (j = 0; j < nbus; ++j) {
			if (getsysinfo(&v, getnode(i), getnode(j), B)) {
				mtx_b->val[i][j] = v;
				fprintf(fp, "%10.5f ", v);
			} else {
				mtx_b->val[i][j] = 0;
				fprintf(fp, "%10s ", " ");
			}
		}
		fputc('\n', fp);
	}
	fclose(fp);
	fmtxprt(logfp, "B", mtx_b);
	return mtx_b;
}

Mtx minq_cnst(void) {
		Mtx b;
		Mtx invb;

		if ((b = sysmatrix()) == NULL) 
			progend(3, "can't get the matrix\n");
		if ((invb = sinv(b)) == NULL)
			progend(3, "matrix inverse error\n");
		fmtxprt(logfp, "inverse of B", invb);
		mtxfree(b);
		return invb;
}

void lpmsg(lprec *lp, void *userhandle, int msg) {
	fprintf(stderr, "lp msg func called\n");
}

/* minq_makelp - 利用输入矩阵cnst（X 阵）建立LP模型 */
void minq_makelp(Mtx cnst) {
		int i;

		if ((lpsolve = open_lpsolve_lib(LP_LIB)) == NULL)
			progend(3, "can't open shared library %s\n", LP_LIB);
		if (!init_lpsolve(lpsolve))
			progend(3, "can't init the lpsolve library\n");

		if ((lp = make_lp(0, dim)) == NULL)/* 建立 LP 模型 */
			progend(3, "can't create new LP model\n");
		set_add_rowmode(lp, TRUE);
		for (i=1; i <= dim; i++)
			cost[i] = 1.0;
		if (!set_obj_fnex(lp, dim, cost, NULL))	/* 设置线性模型目标函数 */
			progend(3, "can't set object function\n");
		fmtxprt(stderr, "constraint:", cnst);
		for (i=0; i < dim; i++)
			fprintf(stderr, " %d", busno[i]);
		fprintf(stderr, "\n");
		for (i=0; i < dim; i++) {
			set_upbo(lp, i+1, qmax[i]);
			if (!selectrow(cnst, busno[i], &cm[1], dim, busno))
				progend(1, "selectrow: %d busno index wrong\n", busno[i]);
			if (!add_constraint(lp, cm, GE, 0))	/* 添加线性模型约束 */
				progend(1, "add_constraint: can't add constraint\n");
		}
		set_add_rowmode(lp, FALSE);
		put_msgfunc(lp, lpmsg, NULL, MSG_LPFEASIBLE | MSG_MILPFEASIBLE | MSG_MILPBETTER);
}

double voltdiff(void) {
	int i;

	maxdiff = 0;
	for (i = 0; i < dim; i++) {
		vdiff[i] = vmin[i] - vac[i];	/* 计算电压偏差 */
		vcormin[i] += alp * vdiff[i];	/* 初始化/修正电压补偿量 */
		vcormax[i] = vmax[i] - vac[i];	/* 电压补偿量的上限 */
		if (maxdiff < ABS(vdiff[i]))    
			maxdiff = ABS(vdiff[i]);
	}
	return maxdiff;
}

int minq_solve(double eps) {
	int i;
	double *s, rhs;

	if (voltdiff() <= eps)
		return 1;			/* no linear programming needed */
	for (i = 0; i < dim; ++i) {						
		if (vcormin[i] > 0)
			rhs = vcormin[i];
		else {
			writelog("voltage of bus[%d] over limit\n", busno[i]+1);
			rhs = 0;
		}
		set_rh(lp, i+1, rhs);	/* 修改线性模型的RHS */
	}
	write_lp(lp, "minq.lp");          /* 输出线性模型 */

	set_timeout(lp, 1);  /* in 1 second */
	set_verbose(lp, 2);
	switch (solve(lp)) {
	case 0:
		get_ptr_primal_solution(lp, &s);
		totq = *s;
		arraycpy(vc, s+1, dim);			/* 约束值 */
		arraycpy(qgen, s + 1 + dim, dim);	/* 电容补偿量 */		
		break;
	case 1:
		progend(3,"sub-optimal\n");
		break;
	case 2:
		progend(3, "infeasible\n");
		break;
	case 3:
		progend(3, "unbounded\n");
		break;
	case 7:
		progend(3, "time out\n");
		break;
	default:
		progend(3, "other errors\n");
		break;
	}
	return 0;	/* linear programming continues */
}

void minq_log(void) {
	int i;

	writelog("___________________\n");
	writelog("%s\t", "bus no.");		/* 控制节点编号 */
	for (i = 0; i < dim; ++i)
		writelog("%10d ", busno[i]+1);
	writelog("\n");
	writelog("rhs\t");  /* 输出电压补偿量(rhs) */
	writearray(logfp, vcormin, dim);
	writelog("Vac\t");  /* 输出潮流计算电压 */
	writearray(logfp, vac, dim);
	writelog("Vmin-Vac\v"); /* 输出电压偏差 */
	writearray(logfp, vdiff, dim);
	writelog("%s\t%10.5f\n", "max", maxdiff);  /* 输出电压偏差 */
	writelog("kvar comp\v"); /* 输出节点电容补偿量(解) */
	writearray(logfp, qgen, dim);
	writelog("%s\t%10.5f\n", "total comp", totq);
	writelog("linerized volt comp\v");	/* 输出节点电压补偿量(解) */
	writearray(logfp, vc, dim);
	writelog("\n");
}

/* 输出错误信息erroinfo，关闭LOG文件logfp并终止程序，返回rc */
void progend(int rc, char *errorinfo, ...) {
	va_list ap;
	char msgbuf[1024];

	va_start(ap, errorinfo);
	vsprintf(msgbuf, errorinfo, ap);
	va_end(ap);
	if (logfp != NULL) {
		fprintf(logfp, "%s", msgbuf);
	}
	fprintf(stderr, "%s", msgbuf);

	/* minq_end(); */
	exit(rc);
}

void writelog(char *title, ...) {
	va_list ap;

	if (logfp != NULL) {
		va_start(ap, title);
		vfprintf(logfp, title, ap);
		va_end(ap);
	}	
}

void writearray(FILE *fp, double *v, int lim) {
	while (lim-- > 0)
		fprintf(fp, "%10.4f ", *v++);
	fprintf(fp, "\n");
}

void minq_end(void) {
	if (logfp != NULL)
		fclose(logfp);
	if (setfp != NULL)
		fclose(setfp);
	if (lffp != NULL)
		closecdf(lffp);
	if (lp != NULL)
		delete_lp(lp);
	if (lpsolve != NULL)
		close_lpsolve_lib(lpsolve);
}

int strtrim(char *s) {
	char *ps;

	for (ps = s+strlen(s)-1; ps >= s; ps--)
		if (!isspace(*ps))
			break;
	*++ps = '\0';
	return ps - s;
}


int readline(FILE *fp, char *line, int len) {
	char *p;

	while ((p = fgets(line, len, fp)) != NULL)
		if (p[0] != '#' && trim(p))
			break;
	return (p == NULL) ? 0 : 1;
}

#define MAXLINE 100

int readdata(FILE *fp) {
	double data;
	int nodeno;
	char name[50];
	char line[MAXLINE];

	readline(fp, line, MAXLINE);
	sscanf(line, "dimension %d", &dim);
	readline(fp, line, MAXLINE);
	sscanf(line, "alpha %lf", &alp);
	for (nodeno = 0; readline(fp, line, MAXLINE); nodeno++) {
		if (!strcmp(line, "node")) {
			while (readline(fp, line, MAXLINE)) {
				sscanf(line, "%s %lf", name, &data);
				if (!strcmp(name, "no")) {
					/* pwf busno's start from 0 */
					busno[nodeno] = (int) data - 1;
					fprintf(stderr, "busno: %d\n", busno[nodeno]);
				} else if (!strcmp(name, "vmin"))
					vmin[nodeno] = data;
				else if (!strcmp(name, "vmax"))
					vmax[nodeno] = data;
				else if (!strcmp(name, "qmax"))
					qmax[nodeno] = data;
				else if (!strcmp(name, "theta"))
					theta[nodeno] = data;
				else if (!strcmp(name, "end"))
					break;
			}
		} else
			return 0;
	}
	return 1;
}
