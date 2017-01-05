#ifndef H_PWF
#define H_PWF

#include <stdio.h>
#include <math.h>
#include "simplex/mtx.h"

#define PI 3.141592654 
#define MAXNODE 500
#define MAXBRANCH 1000
#define MAXJAC 1000
#define fixzero(x) (((x) <= 1e-8 && (x) >= -1e-8) ? 0. : (x))
#define reducerad(x) (((x) >= 2*PI) ? fmod(x, 2*PI) : (x))

enum { TITLE, BUS, BRANCH, END };
enum nodetype { PQ, MVARPQ, PV, SLACK };
enum branchtype { AC, FT, VT, VTMVAR, VPA }; 
enum nodeflag { 
	VOLTOVER = 	   01,
	VOLTUNDER = 	   02,
	PVTOPQ = 	01000 
};
enum { DP, DQ, DV, DA };

struct comp {
	double x, y;
};

struct jacelm {
	double h, m, n, l;
};

struct node {
	int no;
	int flag;
	int type;
	int nconnect;
	char *name;
	double basekv;
	double loadmw, loadmvar;
	struct comp pw;		/* net power injection p.u. == Pg-Pd + j(Qg-Qd)*/
	struct comp pw_act;	/* actual calculated power */
	struct comp volt;	/* initial nodal voltage p.u. */
	struct comp adm_sh;	/* shunt G, B in p.u. */ 
	struct nodechain *nbr;
	struct comp adm_self;	/* nodal self amditance */ 
	double q_min, q_max;	/* reactive generation limits for PV node p.u. */ 
	double volt_ctl;	/* control voltage for PV node p.u. */
	struct jacelm jelm;		/* Jacobian element */
};

struct branch {
	int no;
	int type;
	struct node *inode;	/* tap side node ptr (tap bus) */
	struct node *jnode;	/* z side  node ptr (z bus) */
	double linechar;	/* total line charging, ref to paper 04075293 */
	struct comp k;		/* transformer final turns ratio (complex) */ 
	struct comp adm_line;	/* normal (line) adm of a transformer */
	struct comp adm_se;	/* branch series admitance of pi circuit */
	struct comp iadm_sh;	/* tap side shunt adm of pi circuit */
	struct comp jadm_sh;	/* z side shunt adm of pi circuit*/
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

extern double basemva;	/* system base */

struct comp makecomp(double, double);
struct comp compadd(struct comp, struct comp);
struct comp compmul(struct comp, struct comp);
struct comp compmns(struct comp, struct comp);
struct comp compdiv(struct comp, struct comp);
struct comp compinv(struct comp);
struct comp comprec(struct comp);
struct comp compcnj(struct comp);
double compscale(struct comp);
double angle(struct comp);

void ycalc(void);
void jacalc(struct node *);
struct comp node_pw(struct node *);
struct comp node_flow(struct node *pn);
struct comp b_flow(int, struct branch *);
int gauss(double *errf, double *errx);
int newt(double *errf, double *errx);
struct nodechain *chainalloc(void);
void chainfree(struct nodechain* p);
struct nodechain *addchain(struct nodechain *pc, struct node *pn, struct branch *pb);
struct nodechain *invchain(struct nodechain *h);
struct node *findnode(int no);
struct nodechain *findnbr(struct node *, struct node *);
void getsize(int*, int*);
struct node *getnode(int);
struct branch *getbranch(int);
struct node *addnode(void);
struct branch *addbranch(void);
void clear(void);
void printjac(void);
void printnode(void);
void printybus(void);
void printlinef(void);
int pf(int lim, double tol, char *method);
void setnodeinfo(int, Elm, char *);
Elm getnodeinfo(int, char *);
Elm getsysinfo(struct node*, struct node*, int name);
void reorder(char *);

int trim(char *s);
int titlescan(char *line);
int titlewrite(FILE *fp);
int busscan(char *line);
int branchscan(char *line);
int writebus(struct node *t, FILE *fp);
int writebranch(struct branch *b, FILE *fp);
struct node *makenode(struct node *);
struct branch* makebranch(struct branch *);
void nodecalc(struct node *);
void branchcalc(struct branch *);
void closecdf(FILE *fp);
int readcdf(FILE *fp);
int writecdf(FILE *fp);

#endif
