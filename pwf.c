#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "pwf.h"
#include "mtx.h"
#include "msg.h"
#include "comp.h"

struct nodechain *chainalloc(void) {
	return malloc(sizeof(struct nodechain));
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

/* invchain: inverse a nodechain */
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
int nnode;	/* number of system nodes */
int npvnode;	/* number of PV nodes */
int npqnode; 	/* number of PQ nodes */
int nbranch;	/* number of system branches */
struct node *all_node[MAXNODE];		/* system nodes buffer */
struct node *pv_node[MAXNODE];		/* PV nodes buffer */
struct node *pq_node[MAXNODE];		/* PQ nodes buffer */
struct node *sl_node;			/* slack node ptr */
struct branch *all_branch[MAXBRANCH];	/* system branches buffer */

void getsize(int *maxnd, int *maxbr) {
	*maxnd = nnode;
	*maxbr = nbranch;
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
	struct node **pt, *t;
	struct comp pw;

	printf("%10s %4s %20s%20s%20s\n", 
		"bus no.", "type", 
		"volt (e + jf)", "power (P + jQ)", "power(计算值)");
	loopnode(pt) {
		t = *pt;
		pw = node_pw(t);
		printf("%s%9d %4d %10.4f%10.4f" "%10.4f%10.4f"
			"%10.4f%10.4f\n", 
			(t->flag & (OVER_VOLTLIM|UNDER_VOLTLIM)) ? "*" : " ",
			t->no, 
			t->type,
			compscale(t->volt), angle(t->volt) * 180 / 3.14159, 
			t->pw.x, t->pw.y, pw.x, pw.y);
	}
}

void printybus(void) {
	struct node **t;
	struct nodechain *pc;

	loopnode(t) { 
		printf("%4d: (%.4f %.4f)\n", 
			(*t)->no, (*t)->adm_self.x, (*t)->adm_self.y);
		for (pc = (*t)->nbr; pc != NULL; pc = pc ->next)
			printf("---+%4d: (%.4f %.4f)\n", pc->n->no, 
				- pc->b->adm_se.x, - pc->b->adm_se.y);
	}
}

void printlinef(void) {
	struct node **t;
	struct nodechain *pc;
	struct comp pw;

	loopnode(t) {
		pw = line_flow(*t);
		printf("%4d: ( %.4f, %.4f)\n", (*t)->no, pw.x, pw.y);
		for (pc = (*t)->nbr; pc != NULL; pc = pc ->next)
			printf("--->%4d: (%.4f %.4f)\n", pc->n->no, 
				pc->pw_f.x, pc->pw_f.y);
	}
}

/* ycalc: assembly ybus data */
void ycalc(void) {
	struct node **t;
	struct branch **p;

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

/* return the nodal active power injection P */
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

/* return nodal reactive power injection Q */
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


/* line_flow: calc each power flow of line terminated at the node pn, 
	return the total power flow */
struct comp line_flow(struct node *pn) {
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

struct node *findnode(int no) {
	static struct node **t = NULL;

	if (t != NULL && (*t)->no == no)
		return *t;
	loopnode(t) {
		if ((*t)->no == no)
			return *t;
	}
	return NULL;
}

struct nodechain *findnbr(struct node *base, struct node *aim) {
	struct nodechain *ch;

	for (ch = base->nbr; ch != NULL; ch = ch->next)
		if (ch->n->no == aim->no)
			return ch;
	return NULL;
}

int checknode(void) {
	double q_gen;
	double vt;
	Size i;
	int outcntl = 0;
	struct node *t;

	for (i = 0; (t=getnode(i)) != NULL; i++) {
		switch (t->type) {
		case PV:
			q_gen = t->loadmvar/basemva + t->pw.y;
			/* out control range */
			if (q_gen < t->q_min) {
				q_gen = t->q_min;
				t->flag |= MIN_QGEN;
			} else if (q_gen > t->q_max) {
				q_gen = t->q_max;
				t->flag |= MAX_QGEN;
			} else {	/* within limit */
				continue;
			}
			t->flag |= PVTOPQ; 	/* PV to PQ */
			t->type = PQ;
			t->pw.y = q_gen - t->loadmvar/basemva; /* fixed Q */
			t->volt = makecomp(1.0, .0);	/* flat start */
			t->vt_min = .95;
			t->vt_max = 1.1;
			outcntl++;
			break;
		case PQ:
			vt = compscale(t->volt);
			if (vt > t->vt_max)
				t->flag |= OVER_VOLTLIM;
			else if (vt < t->vt_min)
				t->flag |= UNDER_VOLTLIM;
			break;
		}
	}
	return outcntl;
}

/* gs: perform one cycle of Gauss-Seidel power flow iteration, 
	set max voltage difference of this iteration, *errf, which
	is max |V(r+1) - V(r)| . 
	*errx is used to unify with newt, *errx  == *errf.
	if error occured, return 0,  on success, return 1 */
int gs(double *errf, double *errx) {
	struct comp v_temp, p_temp;
	struct node *t;
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
			t->pw.y = p_temp.y;
			v_temp = node_volt(t);
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
static struct jacidx pwindex[MAXJAC];

void addpwindex(int i, struct node *t, int ltype, int rtype) {
	if (i < 0 || i >= MAXJAC)
		return;
	pwindex[i].n = t;
	pwindex[i].ltype = ltype;
	pwindex[i].rtype = rtype;
}

/* calculate the Jacobian elements related to node t */
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

/* makeindex - create the system Jacobian's index, 
 put the mismach function's left side in arg,
 and set the max error errf.
 the index has the format:

 left  jacobian  right
 side  elements  side
 ___   ___  ___  ___
 DP    H    M    DA
 DQ    N    L    DV 
 ___   ___  ___  ___
*/
int makeindex(Elm *errf, Elm *arg) {
	struct node *t;
	struct jacidx *jacp;
	Size i, j;

	jacp = pwindex;
	j = 0;
	*errf = 0;
	for (i=0; (t=getnode(i)) != NULL; ++i) {
		if (t->type == SLACK)
			continue;	/* skip slack bus */
		jacalc(t);
		/* both PQ and PV buses have DP */
		addpwindex(j, t, DP, DA);
		arg[j] = t->pw.x - t->pw_act.x;
		j++;
		if (t->type == PQ || t->flag & PVTOPQ) { 
		/* PQ or PV-to-PQ buses have 2 indices: DP, DQ */
			addpwindex(j, t, DQ, DV);
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

void recjacalc(struct node *t) {
	Elm a, b, c, d;
	struct comp yik, vk, yii, vi;
	struct nodechain *ch;

	t->pw_act = node_pw(t);
	yii = t->adm_self;
	vi = t->volt;
	a = b = c = d = 0;
	for (ch=t->nbr; ch != NULL; ch=ch->next) {
		yik = compinv(ch->b->adm_se);
		vk = t->volt;
		a += yik.x * vk.x;
		b += yik.y * vk.y;
		c += yik.x * vk.y;
		d += yik.y * vk.x;
		ch->jelm.j1 = -(yik.x * vi.x + yik.y * vi.y);
		ch->jelm.j2 = yik.y * vi.x - yik.x * vi.y;
		switch (t->type) {
		case PQ:
			ch->jelm.j3 = - ch->jelm.j2;
			ch->jelm.j4 = - ch->jelm.j1;
			break;
		case PV:
			ch->jelm.j5 = ch->jelm.j6 = 0;
			break;
		default:
			break;
		}
	}
	t->jelm.j1 = -(a - b + 2 * yii.x * vi.x);
	t->jelm.j2 = -(c + d) - 2 * yii.x * vi.y;
	switch(t->type) {
	case PQ:
		t->jelm.j3 =   c + d  + 2 * yii.y * vi.x;
		t->jelm.j4 = -(a - b) + 2 * yii.y * vi.y;
		break;
	case PV:
		t->jelm.j5 = -2 * vi.x;
		t->jelm.j6 = -2 * vi.y;
		break;
	default:
		break;
	}
}

Size recmakeindex(Elm *errf, Elm *arg) {
	struct node **t;
	struct comp dpw, volt;
	Elm vs;
	Size j;

	j = 0;
	loopnode(t) {
		if ((*t)->type == SLACK)
			continue;
		recjacalc(*t);
		dpw = compmns((*t)->pw, (*t)->pw_act);
		switch((*t)->type) {
		case PQ:
			addpwindex(j, *t, DPR, DVF);
			arg[j] = dpw.x;
			msg(stderr, "pwf: recmakeindex: %d, dP = %f\n", 
				(*t)->no, arg[j]);
			j++;
			addpwindex(j, *t, DQR, DVE);
			arg[j] = dpw.y;
			msg(stderr, "pwf: recmakeindex: %d, dQ = %f\n", 
				(*t)->no, arg[j]);
			j++;
			break;
		case PV:
			addpwindex(j, *t, DPR, DVF);
			arg[j] = dpw.x;
			msg(stderr, "pwf: recmakeindex: %d, dP = %f\n", 
				(*t)->no, arg[j]);
			j++;
			addpwindex(j, *t, DVS, DVE);
			vs = (*t)->volt_ctl;
			volt = (*t)->volt;
			arg[j] = vs * vs 
				- volt.x * volt.x + volt.y * volt.y;
			msg(stderr, "pwf: recmakeindex: %d, dV^2 = %f\n", 
				(*t)->no, arg[j]);
			j++;
			break;
		default:
			break;
		}
	}
	msg(stderr, "pwf: recmakeindex: %4d\n", j);
	norm(arg, j, errf);
	return j;	
}

/* undateindex - update nodal voltages according to the index */
void updateindex(Elm *arg, Size dim, Elm *errx) {
	Size i;
	struct node *t;
	Elm da, dv, ang;
	int rtype;

	norm(arg, dim, errx);
	for (i=0; i < dim; i++) {
		t = pwindex[i].n;
		rtype = pwindex[i].rtype;
		switch(rtype) {
		case DA:
			da = arg[i];
			t->volt = compmul(t->volt, makecomp(cos(da), sin(da)));
			msg(stderr, "pwf: updateindex: da is %g\n", arg[i]);
			break;
		case DV:
			dv = arg[i] + compscale(t->volt);
			ang = angle(t->volt);
			t->volt = makecomp(dv*cos(ang), dv*sin(ang));
			msg(stderr, "pwf: updateindex: dv is %g\n", arg[i]);
			break;
		case DVE:
			t->volt.x -= arg[i];
			msg(stderr, "pwf: updateindex: de is %g\n", arg[i]);
			break;
		case DVF:
			t->volt.y -= arg[i];
			msg(stderr, "pwf: updateindex: df is %g\n", arg[i]);
			break;
		default:
			break;
		}
	}
}

int jactype(int ltype, int rtype) {
	switch (ltype) {
	case DP: 
		if (rtype == DA)
			return H;
		else
			return M;
		break;
	case DQ:
		if (rtype == DA)
			return N;
		else
			return L;
		break;
	case DPR:
		if (rtype == DVE)
			return J1;
		else
			return J2;
		break;
	case DQR:
		if (rtype == DVE)
			return J3;
		else
			return J4;
		break;
	case DVS:
		if (rtype == DVE)
			return J5;
		else
			return J6;
		break;
	default:
		msg(stderr,"pwf: getjac: wrong type in pwindex\n");
		return -1;
	}
}

void makejac(int dim) {
	Size i, j;
	struct jacidx *pi, *pj;

	if (jac != NULL)
		mtxfree(jac);
	jac = mtxalloc(dim, dim);

	for (i = 0; i < dim; i++) {
		pi = &pwindex[i];
		for (j = 0; j < dim; j++) {
			pj = &pwindex[j];
			jac->val[i][j]=getsysinfo(pi->n, pj->n,
				jactype(pi->ltype, pj->rtype));
		}
	}
}

void printjac(void) {
	if (jac != NULL)
		mtxprint("System Jacobian", jac);
}

/* update nodal power */
void updatenp(void) {
	int i;
	struct node *t;
	struct comp pw;

	for (i = 0; (t=getnode(i)) != NULL; ++i)
		switch (t->type) {
		case PV:
			pw = node_pw(t);
			t->pw.y = pw.y;
			msg(stderr, "pwf: update: Qgen of %d is %g\n",
				t->no, pw.y);
			break;
		case SLACK:
			t->pw = node_pw(t);
			break;
		default:
			break;
		}
}

static double df[MAXJAC];  /* mismatch functions' left and right terms */
static int indx[MAXJAC];  /* swap index of ludcmp & lubksb */
static int nrtype;		/* polar form or rectangular form */
enum { POLAR, RECT };

/* nr: Newton-Raphson power flow
	errf - max of dP, dQ
	errx - max of dA, dV or de , df
   return values:
	1 - success
*/
int nr(double *errf, double *errx) {
	Size d;
	Elm det;
	Mtx m;

	d = (nrtype == POLAR)
		? makeindex(errf, df) : recmakeindex(errf, df);
	makejac(d);
	m = mtxdup(jac);
	ludcmp(m, indx, &det);
	lubksb(m, indx, df);
	mtxfree(m);
	updateindex(df, d, errx);
	updatenp();
	return 1;
}

/* pf: power flow calculations
	lim - interation limit
	tol - tolerance
	method - 
		'n' - polar Newton-Raphson
		'r' - rectangular Newton-Raphson
		'g' - Gauss-Seidal
	ischeck -
		1 - regulated load flow
		0 - unregulated load flow
   return values:
	<= lim - convergence
	>  lim - none convergence
	-1 - errors
*/
int pf(int lim, double tol, char *method, int ischeck) {
	int i;
	double errf, errx;
	int (*powerf)(double *, double *);
	
	switch(*method) {
	case 'n':
		powerf = nr;
		nrtype = POLAR;
		break;
	case 'r':
		powerf = nr;
		nrtype = RECT;
		break;
	case 'g':
		powerf = gs;
		break;
	default:
		return -1;
	}

	reorder(all_node, nnode, NODELINES);
	for (i = 1; i <= lim; i++) {
		powerf(&errf, &errx);
		if (ischeck && checknode() > 0)
			continue;
		if (errf <= tol)
			break;
	}
	reorder(all_node, nnode, NODENO);
	return i;
}

Elm getnodeinfo(struct node *t, int name) {
	switch(name) {
		case VOLT:
			return compscale(t->volt);
			break;
		case ANG:
			return angle(t->volt);
			break;
		case QGEN:
			return t->pw.y + t->loadmvar/basemva;
			break;
		default:
			return 0;
	}
	return 0;
}

void setnodeinfo(struct node *t, Elm result, int name) {
	Elm an, volt;

	switch (name) {
	case VOLT:
		an = angle(t->volt);
		t->volt = makecomp(result*cos(an), result*sin(an));
		break;
	case ANG:
		volt = compscale(t->volt);
		t->volt = makecomp(volt*cos(result), volt*sin(result));
		break;
	case QGEN:
		t->pw.y = result - t->loadmvar/basemva;	
		break;
	case QGENINC:
		t->pw.y += result;
		break;
	default:
		break;
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
	case H:
		sen = quad.h;
		break;
	case L:
		sen = quad.l;
		break;
	case M:
		sen = quad.m;
		break;
	case N:
		sen = quad.n;
		break;
	case J1:
		sen = quad.j1;
		break;
	case J2:
		sen = quad.j2;
		break;
	case J3:
		sen = quad.j3;
		break;
	case J4:
		sen = quad.j4;
		break;
	case J5:
		sen = quad.j5;
		break;
	case J6:
		sen = quad.j6;
		break;
	case G:
		sen = cval.x;
		break;
	case B:
		sen = cval.y;
		break;
	default:
		msg(stderr, "pwf: getsysinfo: wrong type %c\n", name);
		return 0;
	}
	return sen;
}

int nodecmp1(const void *t1, const void *t2) {
	return (*(struct node * const *)t1)->nconnect -
		(*(struct node * const *)t2)->nconnect;
}

/* nodecmp2 - compare the nodes according to node no. */
int nodecmp2(const void *t1, const void *t2) {
	return (*(struct node * const *)t1)->no -
		(*(struct node * const *)t2)->no;
}

/* nodecmp3 - compare the nodes according to node types */
int nodecmp3(const void *t1, const void *t2) {
	int t1type, t2type;

	t1type = (*(struct node * const *)t1)->type;
	t2type = (*(struct node * const *)t2)->type;
	return nodetype(t1type) - nodetype(t2type);
}

int nodetype(int ntype) {
	if (ntype == PV) 
		return 0;
	else if (ntype == PQ)
		return 1;
	else if (ntype == SLACK)
		return 2;
	else
		return 3;
}

/* reorder - sort the node list according to type
 nodes:
	all_node or pv_node or pq_node
 type:
	NODELINES  based on connection lines of a node
	NODENO	   based on node no
*/
void reorder(struct node **nodes, int lim, int type) {
	int (*cmpf) (const void *, const void *);

	switch (type) {
	case NODELINES:
		cmpf = nodecmp1;
		break;
	case NODENO:
		cmpf = nodecmp2;
		break;
	case NODETYPE:
		cmpf = nodecmp3;
		break;
	default:
		msg(stderr, "wrong order type\n");
		return;
	}
	qsort(nodes, lim, sizeof(nodes[0]), cmpf);
}

void pvpqsl(void) {
	int i;
	struct node *t;

	npvnode = npqnode = 0;
	for (i = 0; (t=getnode(i)) != NULL; i++)
		switch(t->type) {
		case SLACK:
			sl_node = t;
			break;
		case PV:
			pv_node[npvnode++] = t;
			break;
		case PQ:
			pq_node[npqnode++] = t;
			break;
		default:
			break;
		}
	msg(stderr, "pwf: reorder: pv = %d\n", npvnode);
	msg(stderr, "pwf: reorder: pq = %d\n", npqnode);
	msg(stderr, "pwf: reorder: total = %d\n", nnode);
}

/* flat start of power flow with one cycle of gs iteration */
void flatstart(void) {
	struct node **t;
	double ef, ex, vctl; 
	struct comp volt;

	pvpqsl();
	volt = sl_node->volt;
	loopnode(t) {
		switch((*t)->type) {
		case PQ:
			(*t)->volt = volt; /* set to Slack node's voltage */
			break;
		case PV:
			vctl = (*t)->volt_ctl;
			setnodeinfo(*t, vctl, VOLT);
			break;
		default:
			break;
		}
	}
	gs(&ef, &ex);
}
