#include <stdio.h>
#include <lp_lib.h>
#include "mtx.h"

#define MAXVAR 50
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define LIMT 50
#define TOLR 1e-2

extern FILE *lffp, *logfp, *setfp;

void readvarin(void);
void writevarout(void);
int minq_lf(void);
void minq_start(char *casename, char *logname, char *detaname);
int minq_solvelp(double eps);
Mtx minq_cnst(void);
void minq_makelp(Mtx cnst);
void minq_log(void);
int minq_solve(double eps);
void minq_end(void);
int strtrim(char *);
void writelog(char *title, ...);
void writearray(double *v, int lim);
void progend(int, char *title, ...);
int readline(FILE *fp, char *line, int len);
int readdata(FILE *fp);
int scannode(FILE *fp);
Mtx sysmatrix(void);
