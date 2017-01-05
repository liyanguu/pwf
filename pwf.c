#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include "pwf.h"
#include "simplex/mtx.h"
#include "msg.h"

void print_except(void) {
	int res = fetestexcept(FE_ALL_EXCEPT);

	if (res & FE_INEXACT)
		printf("warning: floating point: loss precision\n");
	if (res & FE_UNDERFLOW)
		printf("warning: floating point: underflow\n");
	if (res & FE_OVERFLOW)
		printf("warning: floating point: overflow\n");
	if (res & FE_DIVBYZERO)
		printf("warning: floating point: div by zero\n");
	if (res & FE_INVALID)
		printf("warning: floating point: invalid operation\n");
	feclearexcept(FE_ALL_EXCEPT);
}

/* makecomp: return a complex number (x, y) */
struct comp makecomp(double x, double y) {
	struct comp temp; 

	temp.x = x;
	temp.y = y;
	return temp;
}

/* compadd: add complex a and complex b, return the sum */
struct comp compadd(struct comp a, struct comp b) {
	a.x += b.x;
	a.y += b.y;
	return a;
}

/* compmns: minus complex a by complex b, return the result */
struct comp compmns(struct comp a, struct comp b) {
	a.x -= b.x;
	a.y -= b.y;
	return a;
}

/* compdiv: divide complex a by complex b, 
	if b = (0, 0) return (0, 0), else return the result */
struct comp compdiv(struct comp a, struct comp b) {
	struct comp temp;
	double buttom;
	
	buttom = b.x * b.x + b.y * b.y;
	if (buttom == 0)
		return makecomp(.0, .0);
	else {
		temp.x = (a.x * b.x + a.y * b.y) / buttom;
		temp.y = (a.y * b.x - a.x * b.y) / buttom;
		return temp;
	}
}

/* compmul: multiply complex a and complex b, return the result */
struct comp compmul(struct comp a, struct comp b) {
	struct comp temp;

	temp.x = a.x * b.x - a.y * b.y;
	temp.y = a.x * b.y + a.y * b.x;
	return temp;
}

/* compmuls: multiply complex a and scaler s, return the result */
struct comp compmuls(struct comp a, double s) {
	a.x *= s;
	a.y *= s;
	return a;
}

/* compinv: return the inverse of complex a = (x, y), i.e. (-x, -y) */
struct comp compinv(struct comp a) {
	a.x = - a.x;
	a.y = - a.y;
	return a;
}

/* compcnj: return the conjunction of complex a = (x, y), i.e. (x, -y) */
struct comp compcnj(struct comp a) {
	a.y = - a.y;
	return a;
}

/* comprec: return the reciprocal of complex a */
struct comp comprec(struct comp a) {
	struct comp temp;
	
	temp = makecomp(1.0, .0);
	return compdiv(temp, a);
}	

/* compscale: return the scale of complex a */
double compscale(struct comp a) {
	double res;

	res = (a.x * a.x + a.y * a.y);
	res = sqrt(res);
	return res;
}

/* angle: return the angle of complex a in rad, 
	from -PI to PI */
double angle(struct comp a) {
	if (a.x >= .0)
		return atan(a.y / a.x);
	else
		return atan(a.y / a.x) + ((a.y >= .0) ? PI : -PI); 
}

/* compdif: return the difference of complex a and complex b,
	i.e. the scale of a - b */ 
double compdif(struct comp a, struct comp b) {
	return compscale(compmns(a, b));
}

struct nodechain *chainalloc(void) {
	return (struct nodechain *) malloc(sizeof(struct nodechain));
}

void chainfree(struct nodechain* p) {
	struct nodechain *q;

	for ( ; p != NULL; p = q) {
		q = p->next;
		free(p);
	}
}

/* addchain: add a member to chain pc, with node pn, branch pb */
struct nodechain *addchain(struct nodechain *pc, struct node *pn, struct branch *pb) {
	struct nodechain *p;

	if ((p = chainalloc()) == NULL)
		return NULL; 
	p->n = pn;
	p->b = pb;
	p->pw_f = makecomp(.0, .0);
	p->next = pc;
	return p;
}

struct nodechain *invchain(struct nodechain *h) {
	struct nodechain *prev, *p;

	if (h == NULL)
		return NULL;
	for (prev = h; (p = h->next) != NULL; prev = p) {
		h->next = p->next;
		p->next = prev;
	}
	return prev;
}


double basemva;

int nnode;
int nbranch;	/* system size */
struct node *all_node[MAXNODE];
struct branch *all_branch[MAXBRANCH];

#define loopnode(t) for (t=all_node; t - all_node < nnode; t++)
#define loopbranch(b) for (b=all_branch; b - all_branch < nbranch; b++)

void getsize(int *maxnd, int *maxbr) {
	*maxnd = nnode;
	*maxbr = nbranch;
}

struct node **getnodelist(int n) {
	if (n < 0 || n >= nnode)
		return NULL;
	return &all_node[n];
}

struct node *getnode(int n) {
	if (n < 0 || n >= nnode)
		return NULL;
	else
		return all_node[n];
}

struct node *addnode(void) {
	struct node *t;

	if (nnode >= MAXNODE)
		return NULL;
	else if ((t=malloc(sizeof(*t))) == NULL)
		return NULL;
	all_node[nnode++] = t;
	return t;
}

void nodefree(struct node *t) {
	if (t != NULL) {
		free(t->name);
		chainfree(t->nbr);
		free(t);
	}
}

void branchfree(struct branch *b) {
	if (b != NULL)
		free(b);
}

struct branch *getbranch(int bn) {
	if (bn < 0 || bn >= nbranch)
		return NULL;
	else
		return all_branch[bn];
}

struct branch *addbranch(void) {
	struct branch *b;

	if (nbranch >= MAXBRANCH)
		return NULL;
	else if ((b=malloc(sizeof(*b))) == NULL)
		return NULL;
	all_branch[nbranch++] = b;
	return b;
}

void clear(void) {
	struct node **t;
	struct branch **b;

	loopnode(t)
		nodefree(*t);
	loopbranch(b)
		branchfree(*b);
	nnode = 0;
	nbranch = 0;
}

int newno(int no) {
	struct node *t;
	int i;

	for (i=0; (t=getnode(i))!= NULL; i++)
		if (t->no == no)
			return i;
	return -1;
}

void printnode(void) {
	int i;
	struct node **pt, *t;
	struct comp pw;

	printf("%10s %4s %20s%20s%20s\n", 
		"bus no.", "type", "volt", "P + jQ", "P + jQ(计算值)");
	loopnode(pt) {
		t = *pt;
		pw = node_pw(t);
		printf("%s%9d %4d %10.4f%10.4f" "%10.4f%10.4f"
			"%10.4f%10.4f\n", 
			(t->flag & (VOLTOVER | VOLTUNDER)) ? "*" : " ",
			t->no, 
			t->type,
			compscale(t->volt), angle(t->volt) * 180 / 3.14159, 
			t->pw.x, t->pw.y, pw.x, pw.y);
	}
}

void printybus(void) {
	struct node **t;
	struct nodechain *pc;

	printf("%4s %20s %20s\n", "bus", "self adm", "mutual adm");
	for (t = all_node ; t < all_node + nnode; t++) {
		printf("%4d: ( %.4f %.4f)", 
			(*t)->no, (*t)->adm_self.x, (*t)->adm_self.y);
		for (pc = (*t)->nbr; pc != NULL; pc = pc ->next)
			printf(" ( %.4f %.4f)", - pc->b->adm_se.x, - pc->b->adm_se.y);
		printf("\n");
	}
}

void printlinef(void) {
	struct node **t;
	struct nodechain *pc;
	struct comp pw;

	printf("%4s %20s %20s\n", "bus", "node flow", "line flow");
	for (t = all_node ; t < all_node + nnode; t++) {
		pw = node_flow(*t);
		printf("%4d: ( %.4f, %.4f)", (*t)->no, pw.x, pw.y);
		for (pc = (*t)->nbr; pc != NULL; pc = pc ->next)
			printf(" ( %.4f %.4f)", pc->pw_f.x, pc->pw_f.y);
		printf("\n");
	}
}

/* ycalc: assembly ybus for all nodes */
void ycalc(void) {
	struct node **t;
	struct branch **p;
	struct comp adm;

	loopbranch(p) {
		(*p)->inode->nbr = addchain((*p)->inode->nbr, 
					(*p)->jnode, *p);	
		(*p)->jnode->nbr = addchain((*p)->jnode->nbr, 
					(*p)->inode, *p);
		branchcalc(*p);
	}
	loopnode(t) {
		nodecalc(*t);
	}
}

/* node_pw: return nodal power injection P + j Q */
struct comp node_pw(struct node *t) {
	struct nodechain *pc;
	struct comp p_temp;

	p_temp = compmul(t->volt, t->adm_self);	/* Vi * Yii */
	/* add by Vk * Yik */
	for (pc = t->nbr; pc != NULL; pc = pc->next)
		p_temp = compadd(p_temp, compmul(pc->n->volt,		
			 compinv(pc->b->adm_se)));	/* Yik = -yik */
	p_temp = compmul(p_temp, compcnj(t->volt));	/* mul by Vi* */
	return compcnj(p_temp);		/* P - jQ => P + jQ*/
}

Elm node_p(struct node *t) {
	Elm vi = compscale(t->volt);
	Elm dan;
	Elm powerp;
	struct nodechain *pc;

	powerp = vi * t->adm_self.x;
	for (pc = t->nbr; pc != NULL; pc = pc->next) {
		dan = angle(t->volt) - angle(pc->n->volt);
		powerp += compscale(pc->n->volt)*(
			-pc->b->adm_se.x * cos(dan) +
			-pc->b->adm_se.y * sin(dan));
	}
	return powerp * vi;
}

Elm node_q(struct node *t) {
	Elm vi = compscale(t->volt);
	Elm dan;
	Elm powerq;
	struct nodechain *pc;

	powerq = -vi * t->adm_self.y;
	for (pc = t->nbr; pc != NULL; pc = pc->next) {
		dan = angle(t->volt) - angle(pc->n->volt);
		powerq += compscale(pc->n->volt)*(
			-pc->b->adm_se.x * sin(dan) +
			pc->b->adm_se.y * cos(dan));
	}
	return powerq * vi;
}

/* node_volt: return nodal voltage V = e + jf */
struct comp node_volt(struct node *t) {
	struct nodechain *pc;
	struct comp v_temp;

	v_temp = compdiv(compcnj(t->pw), compcnj(t->volt));	/* Si* / Vi* */
	for (pc = t->nbr; pc != NULL; pc = pc->next)
		v_temp = compadd(v_temp, 	/* add by - Vk * Yik  = Vk * yik */
			 compmul(pc->b->adm_se, pc->n->volt));
	v_temp = compdiv(v_temp, t->adm_self);	/* div by Yii */
	return v_temp;
}


/* node_flow: calc each power flow of line terminated at the node pn, 
	return the total flow */
struct comp node_flow(struct node *pn) {
	struct comp p_temp = makecomp(.0, .0);
	struct nodechain *pc;

	for (pc = pn->nbr; pc != NULL; pc = pc->next) {
		if (pc->b->inode->no == pn->no)	/* node at primary side */
			pc->pw_f = b_flow(1, pc->b);
		else
			pc->pw_f = b_flow(2, pc->b);	
		p_temp = compadd(pc->pw_f, p_temp);
	}	
	return p_temp;
}

/* b_flow: if dir = 1, return primary side line flow injection, 
	if dir = 2, return secondary side injection of branch pb */
struct comp b_flow(int dir, struct branch *pb) {
	struct comp p_temp;
	struct comp v1, v2, adm;

	if (dir == 1) {
		v1 = pb->inode->volt;
		v2 = pb->jnode->volt;
		adm = pb->iadm_sh;
	} else {
		v1 = pb->jnode->volt;
		v2 = pb->inode->volt;
		adm = pb->jadm_sh;
	}
	p_temp = compmns(v1, v2);	/* (Vi - Vj) */
	p_temp = compmul(p_temp, pb->adm_se);	/* (Vi - Vj) * yij */
	p_temp = compadd(p_temp, 		/* add self part */
		 compmul(v1, adm));		/* Vi * yij0 */
	p_temp = compmul(v1, compcnj(p_temp));	/* mul Vi, Pij = Vi * Iij* */
	return p_temp;
}

struct node 
*findnode(int no) {
	static struct node **t = NULL;

	if (t != NULL && (*t)->no == no)
		return *t;
	loopnode(t)
		if ((*t)->no == no)
			return *t;
	return NULL;
}

struct nodechain 
*findnbr(struct node *base, struct node *aim) {
	struct nodechain *ch;

	for (ch = base->nbr; ch != NULL; ch = ch->next)
		if (ch->n->no == aim->no)
			return ch;
	return NULL;
}

int checknode(void) {
	struct comp pw, v_temp;
	double q_gen, ang;
	Size i;
	int outcntl = 0;
	struct node *t;

	for (i = 0; (t=getnode(i)) != NULL; i++) {
		switch (t->type) {
		case  PV:
			pw = node_pw(t);
			q_gen = t->loadmvar/basemva + pw.y;
			/* out control range */
			if (q_gen < t->q_min) {
				q_gen = t->q_min;
				t->flag |= VOLTOVER;
			} else if (q_gen > t->q_max) {
				q_gen = t->q_max;
				t->flag |= VOLTUNDER;
			} else {	/* within limit */
				t->pw.y = pw.y;
				break;
			}
			t->pw.y = q_gen - t->loadmvar/basemva;
			t->volt = makecomp(1.0, .0);	/* flat start */
			t->type = PQ;		/* PV to PQ */
			t->flag |= PVTOPQ;
			outcntl++;
			break;
		case SLACK:
			t->pw = node_pw(t);
			break;
		}
	}
	return outcntl;
}

/* gauss: perform one iteration of gauss-seidel power flow, 
	set max voltage difference of this iteration, *errf, which
	is max |V(r+1) - V(r)| , *errx is used to unify with newt,
	*errx  == *errf.
	if error occured, return 0,  on success, return 1 */
int gauss(double *errf, double *errx) {
	struct comp v_temp, p_temp;
	struct node *t;
	double ang;		/* volt angle of PV node */
	double q_gen;		/* reactive generation by PV node */
	int i;

	*errf = .0;
	for (i=0; (t=getnode(i)) != NULL; i++) {
		switch (t->type) {
		case SLACK:	/* calc the complex power of slack node */
			t->pw = node_pw(t);
			continue;
		case PQ:	/* calc new node voltage of PQ node */
			v_temp = node_volt(t);
			break;
		case PV:
			p_temp = node_pw(t);	/* calc net power P + jQ */
			q_gen = t->loadmvar / basemva + p_temp.y;
			if (q_gen >= t->q_min && q_gen <= t->q_max) {
				t->pw.y = p_temp.y;
				v_temp = node_volt(t);	/* calc V */ 
				ang = angle(v_temp);	/* new angle and volt */
				v_temp = makecomp(t->volt_ctl * cos(ang),
					  	  t->volt_ctl * sin(ang));
			} else {		/* out control range */
				if (q_gen < t->q_min) {
					q_gen = t->q_min;
					t->flag |= VOLTOVER;
				} else {
					q_gen = t->q_max;
					t->flag |= VOLTUNDER;
				}	
				t->pw.y = q_gen - t->loadmvar/basemva;
				t->volt = makecomp(1.0, .0);	/* flat start */
				t->type = PQ;		/* PV to PQ */
				t->flag |= PVTOPQ;
				v_temp = node_volt(t);		/* new volt */
			}
			break;
		default:
			msg(stderr, 
			"pwf: gauss: wrong node type %d\n", t->type);
			return -1;
		}
		if ((*errf) < compdif(v_temp, t->volt))	/* volt diff */
			*errf = compdif(v_temp, t->volt);
		t->volt = v_temp;	/* update new volt */
	}
	*errx = *errf;
	return 0;
}

Mtx jac;
Spm jacs;
static struct jacidx jacobian_index[MAXJAC];

/* calc Jacobian element node by node */
void jacalc(struct node *t) {
	struct nodechain *ch;
	struct comp pw, vi, vm, yii, yim, dum;
	struct jacelm tmp;
	Elm vv;

	t->pw_act = pw = node_pw(t);
	vi = t->volt;
	yii = t->adm_self;
	vv = compscale(vi);
	vv *= vv;
	tmp.h = -pw.y - vv * yii.y;
	tmp.l = pw.y - vv * yii.y;
	tmp.m = pw.x + vv * yii.x;
	tmp.n = pw.x - vv * yii.x;
	t->jelm = tmp;
	for (ch = t->nbr; ch != NULL; ch = ch->next) {
		vm = ch->n->volt;
		yim = compinv(ch->b->adm_se);
		dum = compmul(vm, yim);
		ch->jelm.h = dum.x * vm.y - dum.y * vm.x;
		ch->jelm.m = dum.x * vm.x + dum.y * vm.y;
		ch->jelm.l = ch->jelm.h;
		ch->jelm.n = -ch->jelm.m;
	}
}

/* create the system Jacobian's index, put the mismach function's
	left side in arg, and set the total error errf */
int makeindex(Elm *errf, Elm *arg) {
	struct node *t;
	struct comp pw;
	struct jacidx *jacp;
	Size i, j;

	jacp = jacobian_index;
	j = 0;
	for (i=0; (t=getnode(i)) != NULL; ++i) {
		if (t->type == SLACK)
			continue;	/* skip slack bus */
		jacalc(t);
		/* PQ and PV bus */
		jacp[j].ltype = DP;
		jacp[j].rtype = DA;
		jacp[j].n = t;
		arg[j] = t->pw.x - t->pw_act.x;
		j++;
		if (t->type == PQ) { /* PQ node has 2 indies, ie. DP, DQ */
			jacp[j].ltype = DQ;
			jacp[j].rtype = DV;
			jacp[j].n = t;
			arg[j] = t->pw.y - t->pw_act.y;
			j++;
		}
	}
	norm(arg, j, errf);
	for (i=0; i < j; i++)
		msg(stderr, "pwf: makeindex: %4d %12.6f\n", 
			jacp[i].n->no, arg[i]);
	msg(stderr, "pwf: makeindex: %4d\n", j);
	return j; 
}

void updateindex(Elm *arg, Size dim, Elm *errx) {
	Size i;
	struct jacidx *jacp = jacobian_index;
	struct node *t;
	Elm an, vl;

	norm(arg, dim, errx);
	for (i=0; i < dim; i++) {
		t = jacp[i].n;
		an = angle(t->volt);
		vl = compscale(t->volt);
		switch(jacp[i].rtype) {
		case DA:
			an = an + arg[i];
			msg(stderr, "pwf: updateindex: da is %g\n", arg[i]);
			break;
		case DV:
			vl += arg[i] * vl;
			msg(stderr, "pwf: updateindex: dv is %g\n", arg[i]);
			break;
		}
		t->volt = makecomp(vl*cos(an), vl*sin(an));
	}
}

int jactype(Size i, Size j) {
	switch (jacobian_index[i].ltype) {
	case DP: 
		if (jacobian_index[j].rtype == DA)
			return 'h';
		else
			return 'm';
		break;
	case DQ:
		if (jacobian_index[j].rtype == DA)
			return 'n';
		else
			return 'l';
		break;
	default:
		msg(stderr,"pwf: getjac: wrong type jacobian_index[%d]\n",i);
		return -1;
	}
}

void makejac(int dim) {
	Size i, j;
	struct node *ti, *tj;

	if (jac != NULL)
		mtxfree(jac);
	jac = mtxalloc(dim, dim);

	for (i = 0; i < dim; i++) {
		ti = jacobian_index[i].n;
		for (j = 0; j < dim; j++) {
			tj = jacobian_index[j].n;
			jac->val[i][j]=getsysinfo(ti, tj, jactype(i, j));
		}
	}
}

void makejacs(int dim) {
	Size i, j;
	struct node *ti, *tj;

	if (jacs != NULL)
		deletsp(jacs);
	jacs = makesp(dim, dim);

	for (i = 0; i < dim; i++) {
		ti = jacobian_index[i].n;
		for (j = 0; j < dim; j++) {
			tj = jacobian_index[j].n;
			add(jacs, i, j, getsysinfo(ti, tj, jactype(i, j)));
		}
	}
}

void printjac(void) {
	if (jac != NULL)
		mtxprint("System Jacobian", jac);
	if (jacs != NULL)
		printsp("Jacobian Sparse", jacs);
}

/* update: update node voltage and angle according to
	mismatch function root arg */
void update(Elm *arg, Size dim, Elm *errx) {
	int i, j;
	struct node *t;
	double vl, an, dan, dvl;
	struct comp pw;

	norm(arg, dim, errx);
	j = 0;
	for (i = 0; (t=getnode(i)) != NULL; ++i)
		switch (t->type) {
		case PQ:
			dan = arg[j++];
			dvl = arg[j++];
			an = angle(t->volt);
			an += dan;
			vl = compscale(t->volt);
			vl += vl * dvl;
			t->volt = makecomp(vl*cos(an), vl*sin(an));
			msg(stderr, "pwf: update: modify %f, %f of %d\n",
				dan, dvl, t->no);
			break;
		case PV:
			dan = arg[j++];
			an = angle(t->volt);
			an += dan;
			vl = compscale(t->volt);
			t->volt = makecomp(vl*cos(an), vl*sin(an));
			pw = node_pw(t);
			t->pw.y = pw.y;
			msg(stderr, "pwf: update: modify %f of %d\n",
				dan, t->no);
			break;
		case SLACK:
			t->pw = node_pw(t);
			break;
		}
}

static double df[MAXJAC];	/* roots of mismatch functions */
static Size indx[MAXJAC];	/* swap index of ludcmp & lubksb */

int newt(double *errf, double *errx) {
	Size d;
	Elm det;

	d = makeindex(errf, df);
	makejac(d);
	ludcmp(jac, indx, &det);
	lubksb(jac, indx, df);
	updateindex(df, d, errx);
	return checknode();
}

int newts(double *errf, double *errx) {
	Size d;

	d = makeindex(errf, df);
	makejacs(d);
	gauss_jordan(jacs, df);
	updateindex(df, d, errx);
	return checknode();
}

int pf(int lim, double tol, char *method) {
	int i;
	double errf, errx;
	int (*powerf)(double *, double *);
	
	switch(*method) {
	case 'n':
		powerf = newt;
		break;
	case 's':
		powerf = newts;
		break;
	case 'g':
		powerf = gauss;
		break;
	default:
		return -1;
	}

	reorder("nconnect");
	for (i = 1; i <= lim; i++) {
		if (powerf(&errf, &errx) > 0)
			continue;
		if (errf <= tol)
			break;
	}
	reorder("busno");
	return i;
}

int getinfo(Elm *v, char *name) {
	struct node **t;
	Elm *pv = v;

	if (strcmp(name, "volt") == 0)
		loopnode(t)
			*v++ = compscale((*t)->volt);
	else
		return -1;
	return v - pv;
}

Elm getnodeinfo(int n, char *name) {
	struct node *t;

	if ((t=getnode(n)) != NULL) {
		if (strcmp(name, "volt")== 0)
			return compscale(t->volt);
		else if (strcmp(name, "qgen")==0)
			return t->pw.y + t->loadmvar/basemva;
	} else
		return .0;
}

void setnodeinfo(int n, Elm result, char *name) {
	struct node *t;
	Elm an;

	if ((t=getnode(n)) != NULL) {
		if (strcmp(name, "volt")==0) {
			an = angle(t->volt);
			t->volt = makecomp(result*cos(an), result*sin(an));
		} else if (strcmp(name, "qgen")==0)
			t->pw.y = result - t->loadmvar/basemva;	
	}
}

Elm getsysinfo(struct node *ti, struct node *tj, int name) {
	struct jacelm quad;
	struct nodechain *ch;
	static Elm sen;
	struct comp cval;

	if (ti->no == tj->no) {
		quad = ti->jelm;
		cval = ti->adm_self;
	} else if ((ch = findnbr(ti, tj)) != NULL) {
		quad = ch->jelm;
		cval = compinv(ch->b->adm_se);
	} else
		return 0;

	switch (name) {
	case 'h':
		sen = quad.h;
		break;
	case 'l':
		sen = quad.l;
		break;
	case 'm':
		sen = quad.m;
		break;
	case 'n':
		sen = quad.n;
		break;
	case 'g':
		sen = cval.x;
		break;
	case 'b':
		sen = cval.y;
		break;
	default:
		msg(stderr, "pwf: getsysinfo: wrong type %c\n", name);
		return 0;
	}
	return sen;
}

int nodecmp(const void *t1, const void *t2) {
	return (*(struct node * const *)t1)->nconnect -
		(*(struct node * const *)t2)->nconnect;
}

int nodenocmp(const void *t1, const void *t2) {
	return (*(struct node * const *)t1)->no -
		(*(struct node * const *)t2)->no;
}

void reorder(char *name) {
	if (strcmp(name, "busno") == 0)
		qsort(all_node, nnode, sizeof(all_node[0]), nodenocmp);
	else if (strcmp(name, "nconnect") == 0)
		qsort(all_node, nnode, sizeof(all_node[0]), nodecmp);
}
