/* pwf.h */

#ifndef H_PWF
#define H_PWF

#include <stdio.h>
#include <math.h>
#include "comp.h"
#include "mtx.h"

#define MAXNODE 500
#define MAXBRANCH 1000
#define MAXJAC 1000
#define fixzero(x) (((x) <= 1e-8 && (x) >= -1e-8) ? 0. : (x))
#define reducerad(x) (((x) >= 2*PI) ? fmod(x, 2*PI) : (x))
#define loopnode(t) for (t=all_node; t - all_node < nnode; t++)
#define loopbranch(b) for (b=all_branch; b - all_branch < nbranch; b++)

enum cdf_type { TITLE, BUS, BRANCH, END, EOD };
/* node_type & branch_type: 
	IEEE Common Format Specification */
enum node_type { PQ = 0, MVARPQ = 1, PV = 2, SLACK = 3 };
enum branch_type { AC, FT, VT, VTMVAR, VPA }; 
enum node_error_flag { 
	VOLTOVER = 	   01,
	VOLTUNDER = 	   02,
	PVTOPQ = 	01000 
};
enum index_type { DP, DQ, DV, DA, DPR, DQR, DVS, DVE, DVF, 
	H, L, M, N , J1, J2, J3 , J4, J5, J6 };
enum node_info_type { G, B, QGEN, QGENINC, VOLT, VOLTINC, ANG, ANGINC };
enum order_type { NODELINES, NODENO, NODETYPE };

extern double basemva;
extern int nnode;	/* number of system nodes */
extern int npvnode;	/* number of PV nodes */
extern int npqnode; 	/* number of PQ nodes */
extern int nbranch;	/* number of system branches */
extern struct node *all_node[MAXNODE];		/* system nodes buffer */
extern struct node *pv_node[MAXNODE];		/* PV nodes buffer */
extern struct node *pq_node[MAXNODE];		/* PQ nodes buffer */
extern struct node *sl_node;			/* slack node ptr */
extern struct branch *all_branch[MAXBRANCH];	/* system branches buffer */

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

void ycalc(void);
void jacalc(struct node *);
struct comp node_pw(struct node *);
struct comp line_flow(struct node *pn);
struct comp b_flow(int, struct branch *);
int gs(double *errf, double *errx);
int nr(double *errf, double *errx);
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
int pf(int lim, double tol, char *method, int ischeck);
void setnodeinfo(struct node *t, Elm, int name);
Elm getnodeinfo(struct node *t, int name);
Elm getsysinfo(struct node *ti, struct node *tj, int name);
int nodetype(int ntype);
void reorder(struct node **nodes, int lim, int type);
void pvpqsl(void);
void flatstart(void);

int trim(char *s);
int titlescan(void);
int titlewrite(void);
int busscan(void);
int branchscan(void);
int writebus(struct node *t);
int writebranch(struct branch *b);
struct node *makenode(struct node *);
struct branch* makebranch(struct branch *);
void nodecalc(struct node *);
void branchcalc(struct branch *);
void closecdf(FILE *fp);
int readcdf(FILE *fp);
int writecdf(FILE *fp);

#endif
