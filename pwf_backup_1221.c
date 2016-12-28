#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pwf.h"
#include "simplex/mtx.h"

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
	return sqrt(a.x * a.x + a.y * a.y);
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

struct node *getnode(int no) {
	if (no < 0 || no >= nnode)
		return NULL;
	else
		return all_node[no];
}

int addnode(struct node *t) {
	if (nnode >= MAXNODE || t == NULL)
		return -1;
	else {
		all_node[nnode++] = t;
		return nnode-1;
	}
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

int addbranch(struct branch *b) {
	if (nbranch >= MAXBRANCH || b == NULL)
		return -1;
	else {
		all_branch[nbranch++] = b;
		return nbranch-1;
	}
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

	printf("%10s %5s %20s %20s\n", "bus no.", "type", "volt", "P + jQ");
	loopnode(pt) {
		t = *pt;
		printf("%s%9d %5d %10.4f%10.4f %10.4f%10.4f\n", 
			(t->flag & (VOLTOVER | VOLTUNDER)) ? "*" : " ",
			t->no, 
			t->type,
			t->volt.x, t->volt.y, 
			t->pw.x, t->pw.y);
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
			printf(" ( %.4f %.4f)", - pc->b->adm_b.x, - pc->b->adm_b.y);
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

/* node_cal: assembly ybus for all nodes */
void node_cal(void) {
	struct branch **p;
	struct comp adm;

	loopbranch(p) {
		adm = compadd((*p)->inode->adm_self, (*p)->adm_b);	/* add branch adm */
		adm = compadd(adm, (*p)->iadm_sh);			/* and shunt adm */
		(*p)->inode->adm_self = adm;			/* to self adm */
		(*p)->inode->nbr = addchain((*p)->inode->nbr, (*p)->jnode, *p);	/* add nbr*/
		adm = compadd((*p)->jnode->adm_self, (*p)->adm_b);	/* the other side */
		adm = compadd(adm, (*p)->jadm_sh);
		(*p)->jnode->adm_self = adm;
		(*p)->jnode->nbr = addchain((*p)->jnode->nbr, (*p)->inode, *p);
	}
}

/* node_pw: return nodal power injection P + j Q */
struct comp node_pw(struct node *t) {
	struct nodechain *pc;
	struct comp p_temp;

	p_temp = compmul(t->volt, t->adm_self);	/* Vi * Yii */
	for (pc = t->nbr; pc != NULL; pc = pc->next)
		p_temp = compadd(p_temp, 
			 compmul(pc->n->volt,		/* add by Vik * Yik */
			 compinv(pc->b->adm_b)));	/* Yik = -yik */
	p_temp = compmul(p_temp, compcnj(t->volt));	/* mul by Vi* */
	return compcnj(p_temp);		/* P - jQ => P + jQ*/
}

/* node_volt: return nodal voltage V = e + jf */
struct comp node_volt(struct node *t) {
	struct nodechain *pc;
	struct comp v_temp;

	v_temp = compdiv(compcnj(t->pw), compcnj(t->volt));	/* Si* / Vi* */
	for (pc = t->nbr; pc != NULL; pc = pc->next)
		v_temp = compadd(v_temp, 	/* add by - Vk * Yik  = Vk * yik */
			 compmul(pc->b->adm_b, pc->n->volt));
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
	p_temp = compmul(p_temp, pb->adm_b);	/* (Vi - Vj) * yij */
	p_temp = compadd(p_temp, 		/* add self part */
		 compmul(v1, adm));		/* Vi * yij0 */
	p_temp = compmul(v1, compcnj(p_temp));	/* mul Vi, Pij = Vi * Iij* */
	return p_temp;
}

struct node 
*findnode(int no) {
	struct node **t;
	int i;

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

/* gauss: perform one iteration of gauss-seidel power flow, 
	set max voltage difference of this iteration, *errf, which
	is max |V(r+1) - V(r)| , *errx is used to unify with newt,
	*errx  == *errf.
	if error occured, return 0,  on success, return 1 */
void gauss(double *errf, double *errx) {
	struct comp v_temp, p_temp;
	struct node **t;
	double ang;		/* volt angle of PV node */
	double q_gen;		/* reactive generation by PV node */

	*errf = .0;
	for (t = all_node; t < all_node + nnode; t++) {
		switch ((*t)->type) {
		case SLACK:	/* calc the complex power of slack node */
			(*t)->pw = node_pw(*t);
			continue;
		case PQ:	/* calc new node voltage of PQ node */
			v_temp = node_volt(*t);
			break;
		case PV:
			p_temp = node_pw(*t);	/* calc net power P + jQ */
			q_gen = (*t)->loadmvar / basemva + p_temp.y;
			if (q_gen >= (*t)->q_min && q_gen <= (*t)->q_max) {
				(*t)->pw.y = p_temp.y;
				v_temp = node_volt(*t);	/* calc V */ 
				ang = angle(v_temp);	/* new angle and volt */
				v_temp = makecomp((*t)->volt_ctl * cos(ang),
					  	  (*t)->volt_ctl * sin(ang));
			} else {		/* out control range */
				if (q_gen < (*t)->q_min) {
					q_gen = (*t)->q_min;
					(*t)->flag |= VOLTOVER;
				} else {
					q_gen = (*t)->q_max;
					(*t)->flag |= VOLTUNDER;
				}	
				(*t)->pw.y = q_gen - (*t)->loadmvar/basemva;
				(*t)->volt = makecomp(1.0, .0);	/* flat start */
				(*t)->type = PQ;		/* PV to PQ */
				(*t)->flag |= PVTOPQ;
				v_temp = node_volt(*t);		/* new volt */
			}
			break;
		default:
			fprintf(stderr, 
			"gauss: wrong node type %d\n", (*t)->type);
			return;
		}
		if (*errf < compdif(v_temp, (*t)->volt))	/* volt diff */
			*errf = compdif(v_temp, (*t)->volt);
		(*t)->volt = v_temp;	/* update new volt */
	}
	*errx = *errf;
}

static Spm jac;
static struct jacidx jacobian_index[MAXJAC];
static Size jacobian_dim;
static double arg[MAXJAC];	/* roots of mismatch functions */
static Size indx[MAXJAC];	/* ! not used : swap index of ludcmp & lubksb */

/* calc Jacobian element node by node */
void jacalc(void) {
	struct node **t;
	struct nodechain *ch;
	struct comp dum, pw, vi, vm, yii, yim;
	struct jacelm tmp;
	double a, vv;

	for (t = all_node; t < all_node+nnode; ++t) {
		pw = node_pw(*t);
		vi = (*t)->volt;
		yii = (*t)->adm_self;
		vv = compscale(vi)*compscale(vi);
		tmp.h = -pw.y - vv * yii.y;
		tmp.l = pw.y - vv * yii.y;
		tmp.m = pw.x + vv * yii.x;
		tmp.n = pw.x - vv * yii.x;
		(*t)->jelm = tmp;
		for (ch = (*t)->nbr; ch != NULL; ch = ch->next) {
			vm = ch->n->volt;
			yim = compinv(ch->b->adm_b);
			a = angle(vm)-angle(vi)+angle(yim);
			dum = compmul(vi, compmul(vm, yim));
			ch->jelm.h = -compscale(dum) * sin(a);
			ch->jelm.m = compscale(dum) * cos(a);
			ch->jelm.l = ch->jelm.h;
			ch->jelm.n = -ch->jelm.m;
		}
	}
}

/* create the system Jacobian's index, put the mismach function's
	left side in arg, and set the total error errf */
void makeindex(double *errf) {
	struct node *t;
	struct comp pw;
	struct jacidx *jacp;
	Size i, j, dim;

	jacp = jacobian_index;
	*errf = .0;
	j = 0;
	for (i=0; (t=getnode(i)) != NULL; ++i) {
		if (t->type == SLACK)
			continue;	/* skip slack bus */
		pw = node_pw(t);
		jacp->ltype = DP;
		jacp->rtype = DA;
		jacp->n = t;
		arg[j] = t->pw.x - pw.x;
		*errf += fabs(arg[j]);
		jacp++;
		j++;
		if (t->type == PQ) {	/* PQ node has 2 indies */
			jacp->ltype = DQ;
			jacp->rtype = DV;
			jacp->n = t;
			arg[j] = t->pw.y - pw.y;
			*errf += fabs(arg[j]);
			jacp++;
			j++;
		}
	}
	jacobian_dim = dim = jacp - jacobian_index;
	for (jacp = jacobian_index; jacp < jacobian_index+dim; ++jacp)
		fprintf(stderr, "makeindex: %4d %12.6f\n", jacp->n->no, arg[i++]);
	fprintf(stderr, "makeindex: %4d\n", dim);
}

Elm getsysinfo(Size i, Size j, int name) {
	struct jacelm quad;
	struct nodechain *ch;
	struct node *ti, *tj;
	static Elm sen;
	struct comp cval;

	if ((ti=getnode(i))==NULL || (tj=getnode(j))==NULL)
		return NULL;
	if (i == j) {
		quad = ti->jelm;
		cval = ti->adm_self;
	} else if ((ch = findnbr(ti, tj)) != NULL) {
		fprintf(stderr, "getsysinfo: find nbr %d, %d\n",
					ti->no, tj->no);
		quad = ch->jelm;
		cval = compinv(ch->b->adm_b);
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
		fprintf(stderr, "getsysinfo: wrong type %c\n", name);
		return 0;
	}
	fprintf(stderr, "getsysinfo: get %f\n", sen);
	return sen;
}
			
int jactype(Size i, Size j) {
	switch (jacobian_index[i].ltype) {
	case DP: 
		if (jacobian_index[j].rtype == DA)
			return 'h';
		else
			retrun 'm';
	case DQ:
		if (jacobian_index[j].rtype == DA)
			return 'n';
		else
			return 'l';
	default:
		fprintf(stderr,"getjac: wrong type jacobian_index[%d]\n",i);
		return 0;
	}
}

void makejac(void) {
	Elm elm;
	Size i, j, newi, newj;

	if (jac == NULL)
		jac = makesp(jacobian_dim, jacobian_dim);
	for (i = 0; i < jacobian_dim; i++) {
		newi = newno(jacobian_index[i].n->no);
		for (j = 0; j < jacobian_dim; j++) {
			newj = newno(jacobian_index[j].n->no);
			elm=getsysinfo(newi, newj, jactype(i, j);
			if (elm != 0)
				add(jac, i, j, elm);
		}
	}
}

void makejac2(void) {
	Size dim, i, j;

	dim = jacobian_dim;
	if (jac == NULL)
		jac = mtxalloc(dim, dim);
	else if (jac->nrow != dim) {
		mtxfree(jac);
		jac = mtxalloc(dim, dim);
	}
	for (i = 0; i < dim; ++i)
		for (j = 0; j < dim; ++j)
			jac->val[i][j] = getjac(i, j);
}

void printjac(void) {
	printsp("System Jacobian", jac);
}

/* update: update node voltage and angle according to
	mismatch function root arg */
void update(double *errx) {
	int i, j;
	struct node *t;
	double vl, an;
	struct comp pw;

	*errx = .0;
	j = 0;
	for (i = 0; (t=getnode(i)) != NULL; ++i)
		switch (t->type) {
		case PQ:
			an = angle(t->volt);
			*errx += fabs(arg[j]);
			an += arg[j];
			j++;
			vl = compscale(t->volt);
			*errx += fabs(arg[j]);
			vl += vl * arg[j];
			j++;
			t->volt = makecomp(vl*cos(an), vl*sin(an));
			break;
		case PV:
			an = angle(t->volt);
			*errx += fabs(arg[j]);
			an += arg[j];
			j++;
			vl = t->volt_ctl;
			pw = node_pw(t);
			t->pw.y = pw.y;
			t->volt = makecomp(vl*cos(an), vl*sin(an));
			break;
		case SLACK:
			t->pw = node_pw(t);
			break;
		}
}

void newt(double *errf, double *errx) {
	jacalc();
	makeindex(errf);
	makejac();
	gauss_jordan(jac, arg);
	update(errx);
}

void jacobian_asub(Elm *x, Elm *y) {
	Size i, j;

	for (i=0; i < jacobian_dim; i++) {
		y[i] = .0;
		for (j=0; j < jacobian_dim; j++)
			y[i] += getjac(i, j) * x[j];
	}
}

void jacobian_atsub(Elm *x, Elm *y) {
	Size i, j;

	for (i=0; i < jacobian_dim; i++) {
		y[i] = .0;
		for (j=0; j < jacobian_dim; j++)
			y[i] += getjac(j, i) * x[j];
	}
}

double mismatch[MAXJAC];

void newtsparse(double *errf, double *errx) {
	int i;
	double rsq;

	jacalc();
	makeindex(errf);
	for (i = 0; i < jacobian_dim; i++) {
		mismatch[i] = arg[i];
		arg[i] = .0;
	}
	sparse(mismatch, jacobian_dim,
		jacobian_asub, jacobian_atsub, 
		arg, &rsq);
	fprintf(stderr, "newtsparse: residual is %10.6f\n", rsq);
	update(errx);
}

int pf(int lim, double tol, char *method) {
	int i;
	double errf, errx;
	void (*powerfunc)(double *, double *);

	if (*method == 'g')
		powerfunc = gauss;
	else if (*method == 'n')
		powerfunc = newt;
	else if (*method == 's')
		powerfunc = newtsparse;
	else
		return -1;

	for (i = 1; i <= lim; i++) {
		(*powerfunc)(&errf, &errx);
		if (errf <= tol || errx <= tol)
			return 1;
	}
	return 0;
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
