/* pwf.h */

#ifndef HEADER_pwf
#define HEADER_pwf

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>	/* for sin, cos */
#include "comp.h"
#include "msg.h"	/* for 信息输出函数 msg() */
#include "mtx.h"	/* for Elm Size & norm */

#define MAXNODE 500
#define MAXBRANCH 1000
#define MAXJAC 1000
#define myabs(x) ((x) < 0 ? -(x) : (x))
#define fixzero(x) (((x) <= 1e-8 && (x) >= -1e-8) ? 0. : (x))
#define reducerad(x) (((x) >= 2*PI) ? fmod(x, 2*PI) : (x))

enum cdf_type { TITLE, BUS, BRANCH, END, EOD };
/* node_type & branch_type: 
	from IEEE Common Format Specification */
enum node_type { PQ = 0, MVARPQ = 1, PV = 2, SLACK = 3, ALL };
enum branch_type { AC, FT, VT, VTMVAR, VPA }; 
enum node_error_flag { 
	OVER_VOLTLIM  =	   01,
	UNDER_VOLTLIM =	   02,
	MIN_QGEN =	010,
	MAX_QGEN =	020,
	PVTOPQ = 	01000 
};
enum index_type { DP, DQ, DV, DA, DPR, DQR, DVS, DVE, DVF, 
	H, L, M, N , J1, J2, J3 , J4, J5, J6 };
enum node_info_type { G, B, PGEN, QGEN, QGENINC, VOLT, VOLTINC, ANG, ANGINC, VOLT_ANG };
enum order_type { NODELINES, NODENO, NODETYPE };
enum powerflow_type { GS, NR_POL, NR_REC };
enum direction { SEND, RECEIV };

extern double basemva;
extern int pftype;	/* power flow type: GS or NR_POL or NR_REC */

struct jacelm {
	double h, m, n, l;
	double j1, j2, j3, j4, j5, j6;
};

struct node {
	int no;
	int flag;
	int type;
	int nconnect;
	char *name;
	double basekv;
	double pload, qload;	/* 有功、无功负荷 (p.u.) */
	struct comp pw;		/* scheduled power p.u. = Pg-Pd + j(Qg-Qd)*/
	struct comp pw_act;	/* actual calculated power p.u. */
	struct comp volt;	/* initial nodal voltage p.u. */
	struct comp adm_sh;	/* shunt G, B in p.u. */ 
	struct comp adm_self;	/* nodal self amditance */ 
	struct nodechain *nbr;
	double q_max, q_min;	/* reactive generation limits for PV node p.u.*/ 
	double vg_max, vg_min;	/* voltage limits for PV nodes p.u. */
	double volt_ctl;	/* control voltage for PV node p.u. */
	double vt_max, vt_min;	/* voltage limits for PQ node p.u. */
	double var_max, var_min; /* var limits for PQ nodes p.u. */
	struct jacelm jelm;		/* Jacobian element */
};

struct branch {
	int no;
	int type;
	struct node *inode;	/* tap side node ptr (tap bus) */
	struct node *jnode;	/* z side  node ptr (z bus) */
	double linechar;	/* total line charging, ref to paper 04075293 */
	double t;		/* transformer final tap ratio  */ 
	struct comp adm_line;	/* normal (line) adm of a transformer */
	struct comp adm_se;	/* branch series admitance of pi circuit */
	struct comp adm_mut;	/* mutual admitance  = -adm_se */
	struct comp adm_ish;	/* tap side shunt adm of pi circuit */
	struct comp adm_jsh;	/* z side shunt adm of pi circuit*/
};

struct nodechain {		/* chain of inter-connected nodes */
	struct nodechain *next;
	struct node *n;
	struct branch *b;
	struct comp pw_f;
	struct jacelm jelm;
};

struct jacidx {
	struct node *n;
	int ltype;
	int rtype;
};

struct nodechain *chainalloc(void);
void nodechainfree(struct nodechain* p);
struct nodechain *addchain(struct nodechain *pc, struct node *pn, struct branch *pb);
struct nodechain *invchain(struct nodechain *h);
struct node *addnode(void);
struct branch *addbranch(void);
struct node *nodealloc(void);
struct branch *branchalloc(void);
void nodefree(struct node *);
void branchfree(struct branch *);

void getsize(int*, int*);
int getnnode(void);
struct node *getnode(int nodeno);
struct branch *getbranch(int branchno);
void loopnode(int);
struct node *nextnode(void);
void loopbranch(void);
struct branch *nextbranch(void);
void pfclear(void);
void pvpqsl(void);
void flatstart(double *ef);

void ycalc(void);
struct comp node_pw(struct node *);
struct comp checknodepw(struct node *);
struct comp line_flow(struct node *);
struct comp b_flow(int, struct branch *);
int gs(double *errf);
struct node *findnode(int no);
struct nodechain *findnbr(struct node *, struct node *);
void printnode(void);
void printybus(void);
void printlinef(void);
void setnodeinfo(struct node *t, int name, ...);
double getnodeinfo(struct node *t, int name);
int getsysinfo(double *val, struct node *ti, struct node *tj, int name);
int jactype(int ltype, int rtype);
int getjac(double *val, int i, int j);
struct jacelm jacalc(struct node *);
struct jacelm jachaincalc(struct node *, struct nodechain *);
void recjacalc(struct node *);
void updatejac(struct node *);
int makeindex(double *ef, double **arg);
int polmakeindex(double *ef, double **arg);
int recmakeindex(double *ef, double **arg);
void updateindex(double *ex);
int checknode(void);

int pf(char *method, int lim, double tol, int ischeck);
void printjac(void);

int trim(char *s);
int titlescan(char *buf);
int busscan(char *buf);
int branchscan(char *buf);
void writetitle(char *buf);
void writebus(char *buf, struct node *pn);
void writebranch(char *buf, struct branch *pb);
struct node *makenode(struct node *);
struct branch* makebranch(struct branch *);
struct comp nodecalc(struct node *);
struct comp branchcalc(struct branch *);
void nodebranch(struct node *, struct branch *);
void closecdf(FILE *fp);
int readcdf(FILE *fp);
int writecdf(FILE *fp);

#endif
