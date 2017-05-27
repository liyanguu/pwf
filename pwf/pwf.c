/* 文件名：pwf.c
 * 描述：辅助潮流计算的基本计算函数
   用法：include "pwf.h"
   更改记录：
   2017-4-19	检查
   2017-5-11	检查
   2017-5-13~14 修正部分函数
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>	/* for sin, cos */
#include "pwf.h"
#include "msg.h"
#include "mtx.h"	/* for Elm Size & norm */

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

	if ((p = chainalloc()) != NULL) {
		p->n = pn;
		p->b = pb;
		p->pw_f = makecomp(.0, .0);
		p->next = pc;
	}
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

static int nodeind;

void loopnode(void) {
    	nodeind = 0;
}

struct node *nextnode(void) {
    	if (nodeind < 0 || nodeind >= nnode)
	    return NULL;
	return all_node[nodeind++];
}

static int branchind;

void loopbranch(void) {
    	branchind = 0;
}

struct branch *nextbranch(void) {
    	if (branchind < 0 || branchind >= nbranch)
	    return NULL;
	return all_branch[branchind++];
}

void getsize(int *maxnd, int *maxbr) {
	*maxnd = nnode;
	*maxbr = nbranch;
}

int getnnode(void) {
	return nnode;
}

struct node *addnode(void) {
	struct node *t;

	if (nnode >= MAXNODE 
	|| (t=nodealloc()) == NULL)
		return NULL;
	return all_node[nnode++] = t;
}

struct node *nodealloc(void) {
	return malloc(sizeof(struct node));
}

void nodefree(struct node *t) {
	if (t != NULL) {
		free(t->name);
		chainfree(t->nbr);
		free(t);
	}
}

struct branch *branchalloc(void) {
	return malloc(sizeof (struct branch));
}

void branchfree(struct branch *b) {
	if (b != NULL)
		free(b);
}

struct branch *addbranch(void) {
	struct branch *b;

	if (nbranch >= MAXBRANCH)
		return NULL;
	else if ((b=branchalloc()) == NULL)
		return NULL;
	all_branch[nbranch++] = b;
	return b;
}

void clear(void) {
	struct node *pn;
	struct branch *pb;

	for (loopnode(); (pn = nextnode()) != NULL; )
		nodefree(pn);
	for (loopbranch(); (pb = nextbranch()) != NULL; )
		branchfree(pb);
	nnode = 0;
	nbranch = 0;
}

int newno(int no) {
	struct node *t;
	int i;

	for (i=0; (t=getnode(i)) != NULL; i++) {
		if (t->no == no) {
			return i;
		}
	}
	return -1;
}

void printnode(void) {
	struct node *pn;
	struct comp pw;
	double pgen, qgen, pload, qload;

	pgen = qgen = pload = qload = 0;

	printf("%10s %4s %20s%20s%20s\n", 
		"bus no.", "type", 
		"volt (e + jf)", "power (P + jQ)", "power(计算值)");
	for (loopnode(); (pn = nextnode()) != NULL; ) {
		pw = pn->pw_act;
		pgen += pn->pw.x + pn->pload;
		qgen += pn->pw.y + pn->qload;
		pload += pn->pload;
		qload += pn->qload;
		printf("%s%9d %4d %10.4f%10.4f" "%10.4f%10.4f" "%10.4f%10.4f\n", 
			(pn->flag & (OVER_VOLTLIM|UNDER_VOLTLIM)) ? "*" : " ",
			pn->no, 
			pn->type,
			compscale(pn->volt), angle(pn->volt) * 180 / 3.14159, 
			pn->pw.x, pn->pw.y, pw.x, pw.y);
	}
	printf("%10s: %.4f\n", "有功发电", pgen);
	printf("%10s: %.4f\n", "有功负荷", pload);
	printf("%10s: %.4f\n", "无功发电", qgen);
	printf("%10s: %.4f\n", "无功负荷", qload);
	printf("%10s: %.4f\n", "有功网损", pgen - pload);
	printf("%10s: %.4f\n", "无功网损", qgen - qload);
}

void printybus(void) {
	struct node *pn;
	struct nodechain *pc;

	for (loopnode(); (pn = nextnode()) != NULL; ) { 
		printf("%4d: (%.4f %.4f)\n", 
			pn->no, pn->adm_self.x, pn->adm_self.y);
		for (pc = pn->nbr; pc != NULL; pc = pc ->next)
			printf("---+%4d: (%.4f %.4f)\n", pc->n->no, 
				- pc->b->adm_se.x, - pc->b->adm_se.y);
	}
}

void printlinef(void) {
	struct node *pn;
	struct nodechain *pc;
	struct comp pw;

	for (loopnode(); (pn = nextnode()) != NULL; ) {
		pw = line_flow(pn);
		printf("%4d: ( %.4f, %.4f)\n", pn->no, pw.x, pw.y);
		for (pc = pn->nbr; pc != NULL; pc = pc ->next)
			printf("--->%4d: (%.4f %.4f)\n", pc->n->no, 
				pc->pw_f.x, pc->pw_f.y);
	}
}

/* ycalc: form the system Ybus */
/* ycalc: 形成系统的Ybus, 即计算节点的自导纳 */
void ycalc(void) {
	struct node *np;

	for (loopnode(); (np = nextnode()) != NULL; nodecalc(np))
	    	;
}

/* node_pw: return calculated nodal power injection Si = P + j Q according to 
	pre-defined node voltages */
struct comp node_pw(struct node *t) {
	struct nodechain *pc;
	struct comp p_temp;

	p_temp = compmul(t->volt, t->adm_self);	/* Vi * Yii */
	/* add by Vm * Yim */
	for (pc = t->nbr; pc != NULL; pc = pc->next)
		p_temp = compadd(p_temp, compmul(pc->n->volt, /* Vm * Yim */		
			 pc->b->adm_mut));
	p_temp = compmul(p_temp, compcnj(t->volt));	/* mul by Vi* */
	return t->pw_act = compcnj(p_temp);		/* P - jQ => P + jQ*/
}

/* checknodepw: set the nodeal actural power ,
   return the power mismatch dS = S(sched) - S(actural) */
struct comp checknodepw(struct node *t) {
    	struct comp pw_mis = makecomp(.0, .0);
	Elm node_q(struct node *);

    	if (t->type == PQ) {
		pw_mis = compmns(t->pw, t->pw_act);
	} else if (t->type == PV) {
		t->pw.y = node_q(t);
		pw_mis.x = t->pw.x - t->pw_act.x;
	} else if (t->type == SLACK) {
		t->pw = t->pw_act;
	}
	return pw_mis;
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
	powerp *= vi;
	return powerp;
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
		powerq += compscale(pc->n->volt) *
		    (-pc->b->adm_se.x * sin(dan) +
			pc->b->adm_se.y * cos(dan));
	}
	powerq *= vi;
	return powerq;
}

/* node_volt: return nodal voltage Vi = e + jf according to 
	it's nodal power Si */
struct comp node_volt(struct node *t) {
	struct nodechain *pc;
	struct comp v_temp;

	v_temp = compcnj(compdiv(t->pw, t->volt));	/* Si* / Vi* */
	for (pc = t->nbr; pc != NULL; pc = pc->next) {
		/* add by -Vk * Yik  = Vk * yik */
		v_temp = compadd(v_temp,
			 compmul(pc->b->adm_se, pc->n->volt));
	}
	v_temp = compdiv(v_temp, t->adm_self);	/* div by Yii */
	return v_temp;
}

/* line_flow: calc power flow for each line terminated at the node pn, 
	return the total power flow */
struct comp line_flow(struct node *pn) {
	struct comp p_temp = makecomp(.0, .0);
	struct nodechain *pc;
 
	for (pc = pn->nbr; pc != NULL; pc = pc->next) {
		if (pc->b->inode->no == pn->no) {
			/* node at sending side */
			pc->pw_f = b_flow(SEND, pc->b);
		} else {
			/* node at receiving side */
			pc->pw_f = b_flow(RECEIV, pc->b);	
		}
		p_temp = compadd(pc->pw_f, p_temp);
	}
	return p_temp;
}

/* b_flow: calculate branch power flow,
if dir = SEND, return line flow injection on the sending side , 
if dir = RECEIV, return line flow injection on the receiving side */
struct comp b_flow(int dir, struct branch *pb) {
	struct comp p_temp;
	struct comp v1, v2, adm;

	if (dir == SEND) {
		v1 = pb->inode->volt;
		v2 = pb->jnode->volt;
		adm = pb->adm_ish;
	} else if (dir == RECEIV) {
		v1 = pb->jnode->volt;
		v2 = pb->inode->volt;
		adm = pb->adm_jsh;
	} else {
	    	msg(stderr, "pwf: b_flow: wrong direction\n");
		return makecomp(.0, .0);
	}
	p_temp = compmns(v1, v2);	/* Vi - Vj */
	p_temp = compmul(p_temp, pb->adm_se);	/* (Vi - Vj) * yij */
	p_temp = compadd(p_temp, 		/* add self part */
		 compmul(v1, adm));		/* Vi * yij0 */
	p_temp = compmul(v1, compcnj(p_temp));	/* mul Vi, Pij = Vi * Iij* */
	return p_temp;
}

struct node *findnode(int no) {
	struct node *pn;

	for (loopnode(); (pn = nextnode()) != NULL; ) {
		if (pn->no == no) {
			return pn;
		}
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

/* checknode: check the node limits.
     for PV node, check its Q generation,
     for PQ node, check its voltage. 
     return the number of out-of-limit nodes */
int checknode(void) {
	double q_gen;
	double vt;
	Size i;
	int outcntl;
	struct node *t;

	outcntl = 0;
	for (i = 0; (t=getnode(i)) != NULL; i++) {
		switch (t->type) {
		case PV:
			q_gen = t->qload + t->pw.y;
			/* out control range */
			if (q_gen < t->q_min) {
				q_gen = t->q_min;
				t->flag |= MIN_QGEN;
			} else if (q_gen > t->q_max) {
				q_gen = t->q_max;
				t->flag |= MAX_QGEN;
			} else {	
			    	/* within limit */
				continue;
			}
			t->flag |= PVTOPQ; 	/* PV to PQ */
			t->type = PQ;
			t->pw.y = q_gen - t->qload; /* make Q fixed */
			t->volt = makecomp(1.0, .0);	/* flat start */
			t->vt_min = .95;
			t->vt_max = 1.1;
			outcntl++;
			break;
		case PQ:
			vt = compscale(t->volt);
			if (vt > t->vt_max) {
				t->flag |= OVER_VOLTLIM;
			} else if (vt < t->vt_min) {
				t->flag |= UNDER_VOLTLIM;
			}
			break;
		}
	}
	return outcntl;
}

/* gs: perform one cycle of Gauss-Seidel power flow iteration, 
	set max voltage difference of this iteration, *errf, which
	is max |V(r+1) - V(r)| , r is the no. of iteration . 
	if error occured, return 0,  on success, return 1 */
int gs(double *errf) {
	struct comp v_temp;
	struct node *t;
	double df, ang;
	int i;

	*errf = .0;
	for (i=0; (t=getnode(i)) != NULL; i++) {
		switch (t->type) {
		case SLACK:	/* calc the complex power of slack node */
			node_pw(t);
			checknodepw(t);
			continue;
		case PQ:	/* calc new node voltage of PQ node */
			v_temp = node_volt(t);
			break;
		case PV:
			node_pw(t);	/* calc node power P + jQ */
			checknodepw(t); /* update the reactive power */
			ang = angle(node_volt(t)); 
			/* update voltage angle */
			v_temp = makecomp(t->volt_ctl * cos(ang), 
				t->volt_ctl * sin(ang));
			break;
		default:
			msg(stderr, 
			"pwf: gs: wrong node type %d\n", t->type);
			return 0;
		}
		df = compdif(v_temp, t->volt); /* |Vi(r+1) - Vi(r)| */
		if (*errf < df)	
			*errf = df;
		t->volt = v_temp;	/* update new volt */
	}
	return 1;
}

static struct jacidx pwindex[MAXJAC];
static double dfbuf[MAXJAC];	/* differences of S(scheduled) & S(calced) */
static int ndf;		/* number of elems in dfbuf */

void addpwindex(int i, struct node *t, int ltype, int rtype) {
	if (i >= 0 && i < MAXJAC) {
		pwindex[i].n = t;
		pwindex[i].ltype = ltype;
		pwindex[i].rtype = rtype;
	}
}

/* getjac: calculate the Jacobian's element J[i, j].
   return 1 if the element is nonzero, 0 if the element is zero */
int getjac(double *jelem, int i, int j) {
	struct jacidx pi, pj;

	pi = pwindex[i];
	pj = pwindex[j];
	return getsysinfo(jelem, pi.n, pj.n, jactype(pi.ltype, pj.rtype));
}

void updatejac(struct node *pn) {
	struct nodechain *pc;

	if (pftype == NR_POL) {
		jacalc(pn);
		for (pc = pn->nbr; pc != NULL; pc = pc->next)
			jachaincalc(pn, pc);
	} else if (pftype == NR_REC) 
	    	recjacalc(pn);
}

/* jaccalc: calculate the self Jacobian elements relates to node t */
struct jacelm jacalc(struct node *t) {
	struct comp pw, vi, yii;
	struct jacelm tmp;
	Elm v2;

	pw = node_pw(t);
	vi = t->volt;
	yii = t->adm_self;
	v2 = compscale(vi);
	v2 *= v2;
	tmp.h = -pw.y - v2 * yii.y;
	tmp.l = pw.y - v2 * yii.y;
	tmp.m = pw.x + v2 * yii.x;
	tmp.n = pw.x - v2 * yii.x;
	return t->jelm = tmp;
}

/* jachaincalc: calculate the mutual J elems relates to node t & its nbrs */
struct jacelm jachaincalc(struct node *t, struct nodechain *tnbr) {
	struct comp vi, vm, yim, dum;

	vi = t->volt;
	vm = tnbr->n->volt;
	yim = compinv(tnbr->b->adm_se); /* Yim = -yim */
	dum = compmul(vm, yim);		/* Vm * Yim  = a + jb*/
	tnbr->jelm.h = dum.x * vi.y - dum.y * vi.x; /* H = a fi - b ei */
	tnbr->jelm.m = dum.x * vi.x + dum.y * vi.y; /* M = a ei + b fi */
	tnbr->jelm.l = tnbr->jelm.h;	/* L = H */
	tnbr->jelm.n = -tnbr->jelm.m;
	return tnbr->jelm;
}

void recjacalc(struct node *t) {
	Elm a, b, c, d;
	struct comp yik, vk, yii, vi;
	struct nodechain *ch;

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
	}
}

void updatedf(double *errf, double **args) {
	int i;
	struct node *t;
	struct comp pw_mis;

	for (i = 0; i < ndf; ) {
		t = pwindex[i].n;
		node_pw(t);
		pw_mis = checknodepw(t);
		switch (t->type) {
		case PQ:
			dfbuf[i] = pw_mis.x;	/* dP = P(sch) - P */
			dfbuf[i+1] = pw_mis.y;	/* dQ = Q(sch) - Q */
			i += 2;
			break;
		case PV:
			dfbuf[i] = pw_mis.x;	/* dP */
			i++;
			break;
		}
	}
	*args = dfbuf;
	norm(dfbuf, ndf, errf);
}

/* makeindex: create the system Jacobian's index pwindex, 
  	put the mismach functions' left side values in dfbuf,
	and set the max error errf.
	return the number of columns of the system Jacobian.

 the system matrix has the format :

 left  jacobian   right
 pwid  matrix     pwid
 side  elements   side
 ___   ___  ___    ___
 DP     H    M     DA
 DQ  =  N    L  *  DV 
 ___   ___  ___    ___
*/

int makeindex(Elm *errf, Elm **arg) {
    	return (pftype == NR_POL) ? polmakeindex(errf, arg) :
	    	recmakeindex(errf, arg);
}

int polmakeindex(Elm *errf, Elm **arg) {
	struct node *t;
	Size i, j;
	struct comp pw_mis;

	j = 0;
	*errf = 0;
	for (i=0; (t=getnode(i)) != NULL; ++i) {
		node_pw(t);
		pw_mis = checknodepw(t);
		updatejac(t);
		if (t->type == SLACK)
			continue;	/* skip slack bus */
		/* both PQ and PV buses have index DP */
		addpwindex(j, t, DP, DA);
		dfbuf[j] = pw_mis.x;
		msg(stderr, "pwf: makeindex: %4d %12.6f\n", t->no, dfbuf[j]);
		j++;
		/* PQ or PV-to-PQ buses have 2 indices: DP, DQ */
		if (t->type == PQ) { 
			addpwindex(j, t, DQ, DV);
			dfbuf[j] = pw_mis.y;
			msg(stderr, "pwf: makeindex: %4d %12.6f\n", t->no, dfbuf[j]);
			j++;
		}
	}
	*arg = dfbuf;
	norm(dfbuf, j, errf);	/* errf = absolute maximum of dfbuf */
	msg(stderr, "pwf: makeindex: %4d, %12.6f\n", j, *errf);
	return ndf = j;
}

Size recmakeindex(Elm *errf, Elm **arg) {
	struct node *pn;
	struct comp pw_mis, volt;
	Elm vs;
	Size j;

	j = 0;
	for (loopnode(); (pn = nextnode()) != NULL; ) {
	    	node_pw(pn);
		pw_mis = checknodepw(pn);
		if (pn->type == SLACK)
			continue;
		recjacalc(pn);
		switch(pn->type) {
		case PQ:
			addpwindex(j, pn, DPR, DVF);
			dfbuf[j] = pw_mis.x;
			msg(stderr, "pwf: recmakeindex: %d, dP = %f\n", 
				pn->no, dfbuf[j]);
			j++;
			addpwindex(j, pn, DQR, DVE);
			dfbuf[j] = pw_mis.y;
			msg(stderr, "pwf: recmakeindex: %d, dQ = %f\n", 
				pn->no, dfbuf[j]);
			j++;
			break;
		case PV:
			addpwindex(j, pn, DPR, DVF);
			dfbuf[j] = pw_mis.x;
			msg(stderr, "pwf: recmakeindex: %d, dP = %f\n", 
				pn->no, dfbuf[j]);
			j++;
			addpwindex(j, pn, DVS, DVE);
			vs = pn->volt_ctl;
			vs *= vs;
			volt = pn->volt;
			dfbuf[j] = vs
				- volt.x * volt.x + volt.y * volt.y;
			msg(stderr, "pwf: recmakeindex: %d, dV^2 = %f\n", 
				pn->no, dfbuf[j]);
			j++;
			break;
		}
	}
	msg(stderr, "pwf: recmakeindex: %4d\n", j);
	norm(dfbuf, j, errf);
	*arg = dfbuf;
	return j;
}

/* undateindex: update nodal voltages according to pwindex & dfbuf.
errx is set to the absolute maximum of dfbuf. */
void updateindex(double *errx) {
	Size i;
	struct node *t;
	Elm da, dv, ang, amp;
	int rtype;

	for (i=0; i < ndf; i++) {
		t = pwindex[i].n;
		rtype = pwindex[i].rtype;
		switch(rtype) {
		case DA:
			da = dfbuf[i];
			t->volt = compmul(t->volt, makecomp(cos(da), sin(da)));
			msg(stderr, "pwf: updateindex: da of %4d is %g\n", 
				t->no, da);
			break;
		case DV:
			amp = compscale(t->volt);
			dv = dfbuf[i] * amp;
			amp += dv;
			ang = angle(t->volt);
			t->volt = makecomp(amp*cos(ang), amp*sin(ang));
			msg(stderr, "pwf: updateindex: dv of %4d is %g\n", 
				t->no, dv);
			break;
		case DVE:
			t->volt.x -= dfbuf[i];
			msg(stderr, "pwf: updateindex: de is %g\n", dfbuf[i]);
			break;
		case DVF:
			t->volt.y -= dfbuf[i];
			msg(stderr, "pwf: updateindex: df is %g\n", dfbuf[i]);
			break;
		}
	}
	norm(dfbuf, ndf, errx);
}

int jactype(int ltype, int rtype) {
	switch (ltype) {
	case DP: 
	    	return (rtype == DA) ? H : M;
	case DQ:
		return (rtype == DA) ? N : L;
	case DPR:
		return (rtype == DVE) ? J1 : J2;
	case DQR:
		return (rtype == DVE) ? J3 : J4;
	case DVS:
		return (rtype == DVE) ? J5 : J6;
	default:
		msg(stderr,"pwf: jactype: wrong type in pwindex\n");
		return -1;
	}
}

/* updatenp: update the nodal power */
void updatenp(void) {
	int i;
	struct node *t;

	for (i = 0; (t=getnode(i)) != NULL; ++i) {
		node_pw(t);
		checknodepw(t);
	}
}

Elm getnodeinfo(struct node *t, int name) {
	switch(name) {
		case VOLT:
			return compscale(t->volt);
		case ANG:
			return angle(t->volt);
		case PGEN:
			return t->pw.x + t->pload;
		case QGEN:
			return t->pw.y + t->qload;
		default:
			return 0;
	}
}

void setnodeinfo(struct node *t, int name, ...) {
	double an;
	double volt;
	double result;
	va_list ap;

	va_start(ap, name);
	switch (name) {
	case VOLT_ANG:
		volt = va_arg(ap, double);
		an = va_arg(ap, double);
		t->volt = makecomp(volt * cos(an), volt * sin(an));
		break;
	case VOLT:
		volt = va_arg(ap, double);
		an = angle(t->volt);
		t->volt = makecomp(volt*cos(an), volt*sin(an));
		break;
	case ANG:
		an = va_arg(ap, double);
		volt = compscale(t->volt);
		t->volt = makecomp(volt*cos(an), volt*sin(an));
		break;
	case QGEN:
		result = va_arg(ap, double);
		t->pw.y = result - t->qload;
		break;
	case QGENINC:
		result = va_arg(ap, double);
		t->pw.y += result;
		break;
	default:
		break;
	}
	va_end(ap);
}

int getsysinfo(double *sen, struct node *ti, struct node *tj, int name) {
	struct jacelm quad;
	struct nodechain *ch;
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
		*sen = quad.h;
		break;
	case L:
		*sen = quad.l;
		break;
	case M:
		*sen = quad.m;
		break;
	case N:
		*sen = quad.n;
		break;
	case J1:
		*sen = quad.j1;
		break;
	case J2:
		*sen = quad.j2;
		break;
	case J3:
		*sen = quad.j3;
		break;
	case J4:
		*sen = quad.j4;
		break;
	case J5:
		*sen = quad.j5;
		break;
	case J6:
		*sen = quad.j6;
		break;
	case G:
		*sen = cval.x;
		break;
	case B:
		*sen = cval.y;
		break;
	default:
		msg(stderr, "pwf: getsysinfo: wrong type %c\n", name);
		return 0;
	}
	return 1;
}

int nodecmp1(const void *t1, const void *t2) {
	return (*(struct node * const *)t1)->nconnect -
		(*(struct node * const *)t2)->nconnect;
}

/* nodecmp2: compare the nodes according to node no. */
int nodecmp2(const void *t1, const void *t2) {
	return (*(struct node * const *)t1)->no -
		(*(struct node * const *)t2)->no;
}

/* nodecmp3: compare the nodes according to node types */
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

/* reorder: sort the nodes according to some type.
 nodes:
	all_node or pv_node or pq_node
 type:
	NODELINES  based on # out lines of a node
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

/* flatstart: flat start of power flow calculation with one cycle of GS iteration */
/* flatstart: 潮流计算的平启动方式，即所有节点的相角等于平衡节点的相角，
 * 所有PQ节点的电压等于平衡节点的电压；所有PV节点的电压等于其设定电压。*/
void flatstart(double *ef) {
	struct node **t;
	struct comp volt;
	double ang;

	pvpqsl();
	volt = sl_node->volt;
	ang = angle(volt);
	looppq(t)	/* set PQ Nodes' voltages to Slack node's voltage */
		(*t)->volt = volt; 
	looppv(t) 	/* set PV nodes' voltages to controled voltage */
		setnodeinfo(*t, VOLT_ANG, (*t)->volt_ctl, ang);
	gs(ef);	/* perform one cycle of G-S flow */
}
