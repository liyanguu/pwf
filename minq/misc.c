/* �ļ�����misc.c
 * ������minq �������Ҫ���㺯��
 * 2017-5-16  ���
 */

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

static FILE *logfp;
static FILE *setfp;
static FILE *lffp;
static hlpsolve lpsolve;	/* ptr to the lpsolve shared library */
#define LP_LIB "lpsolve55.so"	/* name of the lpsolve shared lib */
double vmin[MAXVAR];
double vmax[MAXVAR];
double qmax[MAXVAR];	/* qgen limits */
double theta[MAXVAR];	/* power factor limits */
double alp;		/* iteration factor alpha */

REAL vcormin[MAXVAR];
REAL vcormax[MAXVAR];
REAL cm[MAXVAR];
REAL cost[MAXVAR];		/* LP ���룺����ϵ�� */
REAL vc[MAXVAR];        /* LP �����Ԥ�Ƶĵ�ѹ������ */
REAL qgen[MAXVAR];         /* LP ������ڵ��޹��������� */
REAL totq;				/* LP ��������޹����� */
REAL vac[MAXVAR];		/* �������ݣ���������ĵ�ѹֵ */
REAL vdiff[MAXVAR];		/* ���������ѹ����С��ѹ֮ƫ�� */
REAL maxdiff;			/* ����ѹƫ�� */
static int busno[MAXVAR];	/* ʵ�ʿ��ƽڵ��� (psasp ���ڲ����)*/
static int dim;
lprec *lp;

void setalp(double alpha) {
	alp = alpha;
}

void readvarin(void) {
	int i;

	for (i = 0; i < dim; i++) {
		vac[i] = getnodeinfo(getnode(busno[i]), VOLT);
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
		progend(1, "error: minq_start: can't read load flow data %s\n", casename);
	if ((logfp = fopen(logname, "w")) == NULL)  /* �� log �ļ� */
		progend(1, "error: minq_start: can't open log file %s\n", logname);
	if ((setfp = fopen(dataname, "r")) == NULL) /* ���趨�����ļ�*/
		progend(1, "error: minq_start: can't open data file %s\n", dataname);
	if (!readdata(setfp))
		progend(1, "error: minq_start: can't read data file %s\n", dataname);
	if ((lpsolve = open_lpsolve_lib(LP_LIB)) == NULL)
		progend(3, "error: minq_start: can't open shared library %s\n", LP_LIB);
	if (!init_lpsolve(lpsolve))
		progend(3, "error: minq_start: can't init the lpsolve library\n");

	writelog("minq settings:\n");
	writelog("#Buses\t%d\n", dim);
	writelog("alpha\t%10.4f\n", alp);
	writelog("Busno\t");
	for (i = 0; i < dim; ++i)
		writelog("%10d ", busno[i]+1);
	writelog("\n");
	writelog("Vmin\t");
	writearray(vmin, dim);
	writelog("Vmax\t");
	writearray(vmax, dim);
	writelog("Qmax\t");
	writearray(qmax, dim);
	readcdf(lffp);
	minq_lf();
	writelog("start of linear programming\n\n");
}

int minq_lf(void) {
	int it;

	if ((it = pf("n", LIMT, TOLR, 0)) < 0) {
		writelog("load flow error\n");
		return 0;
	}
	if (it > LIMT) {
		writelog("load flow not convergent\n");
		return 0;
	}
	return 1;
}

/* sysmatrix: ������������ȡϵͳ�Ľڵ㵼�ɾ���
  ����ϵͳ�Ľڵ�B�󣬲��������G�� minq.matrix
modify:
   ʹ��Jacobian����Ԫ L[i,j] ����B����м���
*/
Mtx sysmatrix(void) {
	int i, j;
	int nbus;
	double v;
	Mtx mtx_b;
	FILE *fp;

	if ((fp = fopen("minq.matrix", "w")) == NULL) {
		writelog("error: sysmatrix: can't open minq.matrix\n");
		return NULL;
	}
	if ((nbus = getnnode()) == 0) {
		writelog("error: sysmatrix: no load flow data\n");
		return NULL;
	}
	fprintf(fp, "B part of Ybus:\n");
	for (i = 0; i < nbus; ++i) {
		for (j = 0; j < nbus; ++j)
			if (getsysinfo(&v, getnode(i), getnode(j), B))
				fprintf(fp, "%10.5f ", v);
			else
				fprintf(fp, "%10.5f ", .0);
		fputc('\n', fp);
	}

	if ((mtx_b = mtxalloc(nbus, nbus)) == NULL) {
		writelog("error: sysmatrix: can't alloc matrix\n");
		return NULL;
	}
	fprintf(fp, "\nL in Jacobian:\n");
	for (i = 0; i < nbus; ++i) {
		for (j = 0; j < nbus; ++j) {
			if (getsysinfo(&v, getnode(i), getnode(j), L)) {
				mtx_b->val[i][j] = v;
				fprintf(fp, "%10.5f ", v);
			} else {
				mtx_b->val[i][j] = 0;
				fprintf(fp, "%10.5f ", .0);
			}
		}
		fputc('\n', fp);
	}
	fclose(fp);
	return mtx_b;
}

Mtx minq_cnst(void) {
	Mtx b;
	Mtx invb;

	if ((b = sysmatrix()) == NULL) 
		progend(3, "error: minq_cnst: can't get the matrix\n");
	fmtxprt(logfp, "L", b);
	invb = mtxinv(b);
	if (invb == NULL)
		progend(3, "error: minq_cnst: matrix singular\n");
	fmtxprt(logfp, "inverse", invb);
	mtxfree(b);
	return invb;
}

void lpmsg(lprec *lp, void *userhandle, int msg) {
	writelog("lpmsg: lp msg func called\n");
}

/* minq_makelp - �����������cnst��X ��� L-1 �󣩽���LPģ�� */
void minq_makelp(Mtx cnst) {
	int i;

	if ((lp = make_lp(0, dim)) == NULL)/* ���� LP ģ�� */
		progend(3, "error: minq_makelp: can't create new LP model\n");
	set_add_rowmode(lp, TRUE);
	for (i=1; i <= dim; i++)
		cost[i] = 1.0;
	if (!set_obj_fnex(lp, dim, cost, NULL))	/* ��������ģ��Ŀ�꺯�� */
		progend(3, "error: minq_makelp: can't set object function\n");
	for (i=0; i < dim; i++) {
	    	set_lowbo(lp, i+1, 0.01); /* �ýڵ���С�޹�����Ϊ 0.01 p.u */
		set_upbo(lp, i+1, qmax[i]);
		if (!selectrow(cnst, busno[i], &cm[1], dim, busno))
			progend(1, "error: minq_makelp: %d busno index wrong\n", busno[i]);
		if (!add_constraint(lp, cm, LE, 0))	/* �������ģ��Լ�� */
			progend(1, "error: minq_makelp: can't add %dth constraint\n", i);
	}
	set_add_rowmode(lp, FALSE);
	set_timeout(lp, 1);  /* in 1 second */
	set_verbose(lp, 2);
	put_msgfunc(lp, lpmsg, NULL, 
		MSG_PRESOLVE | MSG_LPOPTIMAL | MSG_LPFEASIBLE );
}

double voltdiff(void) {
	int i;
	double max, rhs;

	max = 0;
	for (i = 0; i < dim; i++) {
		vdiff[i] = vmin[i] - vac[i];	/* �����ѹƫ�� */
		if (max < ABS(vdiff[i]))    
			max = ABS(vdiff[i]);
		vcormin[i] += alp * vdiff[i];	/* ��ʼ��/������ѹ������ */
		vcormax[i] = vmax[i] - vac[i];	/* ��ѹ������������ */
		rhs = vcormin[i];
		if (rhs < 0)
			writelog("warning: voltage of bus[%d] over limit\n", busno[i]+1);
		set_rh(lp, i+1, -rhs);	/* �޸�����ģ�͵�RHS */
	}
	return max;
}

int minq_solve(double eps) {
	double *s;

	maxdiff = voltdiff();
	if (maxdiff <= eps)
		return 1;			/* no linear programming needed */
	write_lp(lp, "minq.lp");          /* �������ģ�� */

	switch (solve(lp)) {
	case 0:
		get_ptr_primal_solution(lp, &s);
		totq = *s;
		arraycpy(vc, s+1, dim);			/* Լ��ֵ */
		arraycpy(qgen, s + 1 + dim, dim);	/* ���ݲ����� */		
		break;
	case 1:
		progend(3,"error: minq_solve: sub-optimal\n");
		break;
	case 2:
		progend(3, "error: minq_solve: infeasible\n");
		break;
	case 3:
		progend(3, "error: minq_solve: unbounded\n");
		break;
	case 7:
		progend(3, "error: minq_solve: time out, no solutions\n");
		break;
	default:
		progend(3, "error: minq_solve: unknown errors\n");
		break;
	}
	return 0;	/* linear programming continues */
}

void minq_log(void) {
	int i;

	writelog("___________________\n");
	writelog("%s\t", "bus no.");		/* ���ƽڵ��� */
	for (i = 0; i < dim; ++i)
		writelog("%10d ", busno[i]+1);
	writelog("\n");
	writelog("Vac\t");  /* ������������ѹ */
	writearray(vac, dim);
	writelog("Vmin-V\t"); /* �����ѹƫ�� */
	writearray(vdiff, dim);
	writelog("%s\t%10.4f\n", "max", maxdiff);  /* �����ѹƫ�� */
	writelog("Vcorr\t");  /* �����ѹ������(rhs) */
	writearray(vcormin, dim);
	writelog("kvar\t"); /* ����ڵ���ݲ�����(��) */
	writearray(qgen, dim);
	writelog("%s\t%10.4f\n", "total", totq);
	writelog("linVC\t");	/* ����ڵ��ѹ������(��) */
	writearray(vc, dim);
	writelog("\n");
}

/* ���������Ϣerroinfo���ر�LOG�ļ�logfp����ֹ���򣬷���rc */
void progend(int rc, char *errorinfo, ...) {
	va_list ap;
	char msgbuf[1024];

	va_start(ap, errorinfo);
	vsprintf(msgbuf, errorinfo, ap);
	va_end(ap);
	writelog("%s", msgbuf);
	fprintf(stderr, "%s", msgbuf);
	minq_end();
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

void writearray(double *v, int lim) {
	while (lim-- > 0)
		writelog("%10.4f ", *v++);
	writelog("\n");
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
		if (p[0] != '#' && strtrim(p))
			break;
	return (p == NULL) ? 0 : 1;
}

#define MAXLINE 100

int readdata(FILE *fp) {
	char line[MAXLINE];

	while (readline(fp, line, MAXLINE)) {
		if (sscanf(line, "dimension %d", &dim) == 1)
			continue;
		else if (sscanf(line, "alpha %lf", &alp) == 1)
			continue;
		else if (!strcmp(line, "node")) {
			if (!scannode(fp))
				return 0;
		} else
			return 0;
	}
	return (dim <= 0 || alp <= 0) ? 0 : 1;
}

int scannode(FILE *fp) {
	char name[50];
	char line[MAXLINE];
	double data;
	static int nodeno = 0;

	while (readline(fp, line, MAXLINE)) {
		sscanf(line, "%s %lf", name, &data);
		if (!strcmp(name, "no")) {
			/* pwf busno's start from 0 */
			busno[nodeno] = (int) data - 1;
		} else if (!strcmp(name, "vmin"))
			vmin[nodeno] = data;
		else if (!strcmp(name, "vmax"))
			vmax[nodeno] = data;
		else if (!strcmp(name, "qmax"))
			qmax[nodeno] = data;
		else if (!strcmp(name, "theta"))
			theta[nodeno] = data;
		else if (!strcmp(name, "end")) {
			nodeno++;
			return 1;
		} else {
			fprintf(stderr, "unkown field name: %s\n", name);
			break;
		}
	}
	fprintf(stderr, "incomplete node section\n");
	return 0;
}
