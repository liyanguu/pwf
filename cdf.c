/* change log:
	2017-3-19	modify bus section data, add ctlmin,max
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "pwf.h"
#include "comp.h"

		/* 6 items */
#define TITLEFMT	" %8c %20c %lf %d %1s %30c"
#define TITLEPFMT	" %-8s %-20s %6.1f %4d %1s %-30s\n"

		/* 18 terms */
#define BUSFMT "%d %12c %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d"

#define BUSPFMT "%4d %12s %2d %2d %2d %6.3f%7.2f%9.1f%10.1f%8.1f%8.1f %7.1f %6.3f%8.1f%8.1f%8.1f%8.1f %4d\n"

		/* 21 terms */
#define BRANCHFMT "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf"

#define BRANCHPFMT "%4d %4d%3d%3d %1d %1d%10.5f%10.5f%10.5f%6d%6d%6d%4d %d %7.4g %7.4g%7.4g%7.4g%7.5g %7.4g%7.4g\n"

int section, lineno;

char senddate[9], sendername[21], season[2], id[31];
int year;

int busno, lfano, lzno, bustype, rcbusno;
double finalvolt, finalang, loadmw, loadmvar;
double genmw, genmvar, basekv, volt_ctl;
double opmax, opmin; 	/* operation parameters' limits */
double shcon, shsus;
double ctlmax, ctlmin;	/* control parameters' limits */
char name[13];	/* a bus name field is 6-17, 12 chars */

int tbusno, zbusno;
int circuit, branchtype;
double br, bx, linechar;  
double mvarate1, mvarate2, mvarate3;
int ctlbusno, ctlside;
double tfinalturn, tfinalang, mintap, maxtap;
double stepsize, minvolt, maxvolt;

#define MAXLINE 150
char line[MAXLINE];	/* a cdf line */

int readcdf(FILE *fp) {
	char *p;
	int trim(char *);

	while ((p = fgets(line, MAXLINE, fp)) != NULL) {
		if (p[0] == '#' || trim(p) == 0)
			continue;
		++lineno;
		if (lineno == 1) {
			if (!titlescan()) {
				fprintf(stderr, 
					"readcdf: missing title card\n");
				return 0;
			}
		} else if(strstr(line, "BUS DATA FOLLOWS") == line)
			section = BUS;
		else if(strstr(line, "BRANCH DATA FOLLOWS") == line)
			section = BRANCH;
		else if (strstr(line, "-999") == line
			|| strstr(line, "-99") == line 
			|| strstr(line, "-9") == line)
			section = END;
		else if (strstr(line, "END OF DATA") == line)
			section = EOD;
		else {
			switch (section) {
			case BUS:
				if (!busscan())
					return 0;
				break;
			case BRANCH:
				if (!branchscan())
					return 0;
				break;
			}
		}
	}
	if (section != EOD) {
		fprintf(stderr, "readcdf: missing END OF DATA\n");
		return 0;
	}
	ycalc();
	return 1;
}

/* writecdf: write load flow result in cdf to fp */
int writecdf(FILE *fp) {
	struct node *t;
	struct branch *b;
	time_t mytime;
	struct tm *loctime;
	char month[5];
	char yearbuf[5];
	int m;
	int i, maxnode, maxbranch;

	getsize(&maxnode, &maxbranch);
	time(&mytime);
	loctime = localtime(&mytime);
	strftime(senddate, sizeof(senddate), "%D", loctime);
	strncpy(sendername, "WHU LAB", sizeof(sendername)-1);
	strftime(yearbuf, 5, "%Y", loctime);
	year = atoi(yearbuf);
	strftime(month, 5, "%m", loctime);
	strcpy(season, ((m=atoi(month)) <= 8) ? "S": "W");
	titlewrite();
	fprintf(fp, "%s", line);

	fprintf(fp, "BUS DATA FOLLOWS %30d ITEMS\n", maxnode);
	for (i = 0; (t=getnode(i)) != NULL; i++) {
		writebus(t);
		fprintf(fp, "%s", line);
	}
	fprintf(fp, "-999\n");

	fprintf(fp, "BRANCH DATA FOLLOWS %30d ITEMS\n", maxbranch);
	for (i = 0; (b=getbranch(i)) != NULL; i++) {
		writebranch(b);
		fprintf(fp, "%s", line);
	}
	fprintf(fp, "-999\n");

	fprintf(fp, "END OF DATA\n");

	return !ferror(fp);
}

void closecdf(FILE *fp) {
	fclose(fp);
	lineno = 0;
}

int titlescan(void) {
	int n;

	n = sscanf(line, TITLEFMT, senddate, sendername, &basemva,
			&year, season, id);
	if (n != 6) {
		fprintf(stderr, 
		"error: can't read title in this line:\n%s\n", line);
		return 0;
	} else {
		senddate[8] = '\0';
		sendername[20] = '\0';
		season[1] = '\0';
		id[30] = '\0';
		return 1;
	}
}

int titlewrite(void) {
	return sprintf(line, TITLEPFMT, senddate, sendername, basemva,
			year, season, id);
}

int busscan(void) {
	struct node *pt;

	if (sscanf(line, BUSFMT, 
			&busno, name, 
			&lfano, &lzno, &bustype,
			&finalvolt, &finalang, 
			&loadmw, &loadmvar,
			&genmw, &genmvar, 
			&basekv, &volt_ctl, 
			&opmax, &opmin, 
			&shcon, &shsus, 
			&rcbusno) != 18) {
		fprintf(stderr, 
		"busscan: can't read bus in this line:\n %s\n", line);	
		return 0;
	}
	name[12] = '\0';
	if ((pt = addnode()) == NULL) {
		fprintf(stderr, "busscan: can't add node\n");
		return 0;
	}
	makenode(pt);
	return 1;
}

int branchscan(void) {
	struct branch *pb;

	if (sscanf(line, BRANCHFMT, 
	&tbusno, &zbusno, 
	&lfano, &lzno,
	&circuit, &branchtype, 
	&br, &bx, &linechar,
	&mvarate1, &mvarate2, &mvarate3,
	&ctlbusno, &ctlside, 
	&tfinalturn, &tfinalang,
	&mintap, &maxtap, 
	&stepsize, 
	&minvolt, &maxvolt) != 21) {
		fprintf(stderr, "branchscan: wrong cdf line:\n%s\n",
			line);
		return 0;
	}
	if ((pb=addbranch()) == NULL) {
		fprintf(stderr, "branchscan: can't add branch\n");
		return 0;
	}
	makebranch(pb);
	return 1;
}

int writebus(struct node *t) {
	double genmw, genmvar;
	double pmax, pmin;
	int bustype;
	int nbytes;
	
	bustype = t->type;
	if (bustype == PQ && t->flag & PVTOPQ)
		bustype = PV;
	if (bustype == PV) {
		pmax = t->q_max * basemva;
		pmin = t->q_min * basemva;
	} else if (bustype == PQ) {
		pmax = t->vt_max;
		pmin = t->vt_min;
	}
	genmw = fixzero(t->pw.x * basemva + t->loadmw);
	genmvar = fixzero(t->pw.y * basemva + t->loadmvar);
	nbytes = sprintf(line, BUSPFMT,
		t->no, t->name, 1, 1, bustype,
		compscale(t->volt), angle(t->volt) * 180.0 / PI,
		t->loadmw, t->loadmvar,
		genmw, genmvar,
		t->basekv, t->volt_ctl,
		pmax, pmin,
		t->adm_sh.x, t->adm_sh.y, 0);
	return nbytes;
}

int writebranch(struct branch *b) {
	struct comp imp;
	int nbytes;

	imp = comprec(b->adm_line);
	nbytes = sprintf(line, BRANCHPFMT,
		b->inode->no, b->jnode->no,
		1, 1, 1, b->type,
		imp.x, imp.y, b->linechar, 
		0, 0, 0, 0, 0, b->t, .0,
		.0, .0, .0, .0, .0);
	return nbytes;
}

int trim(char *s) {
	char *p = s + strlen(s) - 1;

	for ( ; p >= s; p--)
		if (!isspace(*p))
			break;
	*++p = '\0';
	return p - s;
}

/* makenode: make a node from cdf bus data */
struct node *makenode(struct node *t) {
	double ang_rad;

	t->no = busno;
	t->flag = 0;
	t->type = bustype;
	t->name = strdup(name);
	t->basekv = basekv;
	t->nbr = NULL;
	t->loadmw = loadmw;
	t->loadmvar = loadmvar;
	if (t->type == PV) {
	/* for PV nodes, opmin = min MVAr, opmax = max mvar */
		t->q_min = opmin / basemva;
		t->q_max = opmax / basemva;
	} else if (t->type == PQ) {
	/* for PQ nodes, minop = min voltage, maxop = max volt limit */
		t->vt_min = opmin;
		t->vt_max = opmax;
		t->var_min = ctlmin;
		t->var_max = ctlmax;
	}
	t->volt_ctl = volt_ctl;
	t->adm_sh = makecomp(shcon, shsus);
	t->adm_self = makecomp(0, 0);
	ang_rad = finalang / 180.0 * 3.141592654;
	t->volt = makecomp(finalvolt * cos(ang_rad), 
			finalvolt * sin(ang_rad));
	t->pw = makecomp((genmw - loadmw)/basemva,
			(genmvar - loadmvar)/basemva);
	return t;
}

struct branch *makebranch(struct branch *b) {
	struct comp imp_b;

	b->inode = findnode(tbusno);	
	b->jnode = findnode(zbusno);
	if (b->inode == NULL || b->jnode == NULL) {
		fprintf(stderr, 
		"error: line without terminal nodes:\n%s\n", line);
		return NULL;
	}
	b->type = branchtype;
	b->linechar = linechar;
	imp_b = makecomp(br, bx);
	b->adm_line = comprec(imp_b);
	b->t = tfinalturn;
	return b;
}

/* calculate the pi-model of transmission lines */
void branchcalc(struct branch *b) {
	struct comp adm; 
	double t;

	switch (b->type) {
	case AC:	/* ac line */
		b->adm_se = b->adm_line;
		b->adm_ish = b->adm_jsh = makecomp(.0, b->linechar/2.0);
		break;
	case FT:	/* fixed tap transformer */
		adm = b->adm_line;
		t = b->t;
		b->adm_se = compmuls(adm, t); /* Tij * yij */
		b->adm_ish = compmuls(adm, t*t - t); /* (Tij ** 2 - Tij)*yij*/
		b->adm_jsh = compmuls(adm, 1 - t);  /* (1 - Tij) * yij */
		break;
	default:
		fprintf(stderr, "cdf error: wrong branch type, %d\n", b->type);
		break;
	}
}

/* calculate the self-admitance of the node */
void nodecalc(struct node *t) {
	struct nodechain *pc;
	struct comp adm;

	adm = t->adm_sh;
	t->nconnect = 0;
	for (pc = t->nbr; pc != NULL; pc = pc->next) {
		adm = compadd(pc->b->adm_se, adm);
		if (t == pc->b->inode)
			adm = compadd(pc->b->adm_ish, adm);
		else
			adm = compadd(pc->b->adm_jsh, adm);
		t->nconnect++;
	}
	t->adm_self = adm;
}

/* initdata: read from cdf file fp and initialize system data ,
	return 1 on success, 0 otherwise
int initdata(FILE *fp) {
	if (readcdf(line, MAXLINE, fp) != TITLE) {
		fprintf(stderr, "cdf error: missing title card\n");
		return 0;
	}
	titlescan(line);
	if (readcdf(line, MAXLINE, fp) != BUS) {
		fprintf(stderr, "cdf error: missing bus card\n");
		return 0;
	}
	while (readcdf(line, MAXLINE, fp) != END)
		if (!busscan(line)) 
			return 0;
	if (readcdf(line, MAXLINE, fp) != BRANCH) {
		fprintf(stderr, "cdf error: missing branch card\n");
		return 0;
	}
	while (readcdf(line, MAXLINE, fp) != END)
		if (!branchscan(line))
			return 0;
	return 1;
}
*/
