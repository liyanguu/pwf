#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "pwf.h"

		/* 6 items */
#define TITLEFMT	" %8c %20c %lf %4c %1s %30c"
#define TITLEPFMT	" %8s %20s %6.1f %4s %1s %30s\n"

		/* 18 terms */
#define BUSFMT "%d %12c %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d"

#define BUSPFMT "%4d %12s %2d %2d %2d %6.3f%7.2f%9.1f%10.1f%8.1f%8.1f %7.1f %6.3f%8.1f%8.1f%8.1f%8.1f %4d\n"

		/* 21 terms */
#define BRANCHFMT "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf"

#define BRANCHPFMT "%4d %4d%3d%3d %1d %1d%10.5f%10.5f%10.5f%6d%6d%6d%4d %d %7.4g %7.4g%7.4g%7.4g%7.5g %7.4g%7.4g\n"

int section, lineno;

char senddate[10], sendername[22], year[6], season[2], id[32];

int busno, lfano, lzno, bustype, rcbusno;
double finalvolt, finalang, loadmw, loadmvar;
double genmw, genmvar, basekv, volt_ctl;
double maxmvar, minmvar, shcon, shsus;
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
		if (trim(p) == 0)
			continue;
		if (++lineno == 1) {
			section = TITLE;
			if (!titlescan(line)) {
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
			section = EOF;
		else
			switch (section) {
			case BUS:
				if (!busscan(line)) {
					fprintf(stderr, "readcdf: can't read bus data\n");
					return 0;
				}
				break;
			case BRANCH:
				if (!branchscan(line)) {
					fprintf(stderr, "readcdf: can't read branch data\n");
					return 0;
				}
				break;
			}
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
	int m;
	int i, maxnode, maxbranch;

	getsize(&maxnode, &maxbranch);
	time(&mytime);
	loctime = localtime(&mytime);
	strftime(senddate, sizeof(senddate), "%D", loctime);
	strncpy(sendername, "WHU LAB", sizeof(sendername)-1);
	strftime(year, sizeof(year), "%Y", loctime);
	strftime(month, 5, "%m", loctime);
	strcpy(season, ((m=atoi(month)) <= 8) ? "S": "W");
	titlewrite(fp);

	fprintf(fp, "BUS DATA FOLLOWS %30d ITEMS\n", maxnode);
	for (i = 0; (t=getnode(i)) != NULL; i++)
		writebus(t, fp);
	fprintf(fp, "-999\n");

	fprintf(fp, "BRANCH DATA FOLLOWS %30d ITEMS\n", maxbranch);
	for (i = 0; (b=getbranch(i)) != NULL; i++)
		writebranch(b, fp);
	fprintf(fp, "-999\n");

	fprintf(fp, "END OF DATA\n");

	return !ferror(fp);
}

void closecdf(FILE *fp) {
	fclose(fp);
	lineno = 0;
}

int titlescan(char *line) {
	if (sscanf(line, TITLEFMT,
			senddate,
			sendername,
			&basemva,
			year,
			season,
			id) != 6)
		return 0;
	else {
		senddate[9] = '\0';
		sendername[21] = '\0';
		year[5] = '\0';
		season[1] = '\0';
		id[31] = '\0';
		return 1;
	}
}

int titlewrite(FILE *fp) {
	fprintf(fp, TITLEPFMT,
			senddate,
			sendername,
			basemva,
			year, 
			season, 
			id);
	return !ferror(fp);
}

int busscan(char *line) {
	struct node *pt;

	if (sscanf(line, BUSFMT, 
			&busno, name, 
			&lfano, &lzno, &bustype,
			&finalvolt, &finalang, 
			&loadmw, &loadmvar,
			&genmw, &genmvar, 
			&basekv, &volt_ctl, 
			&maxmvar, &minmvar, 
			&shcon, &shsus, 
			&rcbusno) != 18) 
		return 0;
	name[13] = '\0';
	if ((pt = addnode()) == NULL) {
		fprintf(stderr, "busscan: can't make node\n");
		return 0;
	}
	makenode(pt);
	return 1;
}

int branchscan(char *line) {
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
		fprintf(stderr, "branchscan: wrong cdf line:\n\t%s\n",
			line);
		return 0;
	}
	if ((pb=addbranch()) == NULL) {
		fprintf(stderr, "branchscan: can't make branch\n");
		return 0;
	}
	makebranch(pb);
	return 1;
}

int writebus(struct node *t, FILE *fp) {
	double genmw, genmvar;
	int bustype;
	
	bustype = t->type;
	if (bustype == PQ && t->flag & PVTOPQ)
		bustype = PV;
	genmw = t->pw.x * basemva + t->loadmw;
	genmvar = t->pw.y * basemva + t->loadmvar;
	genmw = fixzero(genmw);
	genmvar = fixzero(genmvar);
	fprintf(fp, BUSPFMT,
		t->no, t->name, 1, 1, bustype,
		compscale(t->volt), angle(t->volt) * 180.0 / PI,
		t->loadmw, t->loadmvar,
		genmw, genmvar,
		t->basekv, t->volt_ctl,
		t->q_max * basemva,
		t->q_min * basemva,
		t->adm_sh.x, t->adm_sh.y, 0);
    	return !ferror(fp);
}

int writebranch(struct branch *b, FILE *fp) {
	struct comp imp;

	imp = comprec(b->adm_line);
	fprintf(fp, BRANCHPFMT,
		b->inode->no, b->jnode->no,
		1, 1, 1, b->type,
		imp.x, imp.y, b->linechar, 
		0, 0, 0, 0, 0, b->k.x, b->k.y,
		.0, .0, .0, .0, .0);
	return !ferror(fp);
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
	t->q_min = minmvar / basemva;
	t->q_max = maxmvar / basemva;
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
			"cdf error: line without terminal node, %d\n", lineno);
		return NULL;
	}
	b->type = branchtype;
	b->linechar = linechar;
	imp_b = makecomp(br, bx);
	b->adm_line = comprec(imp_b);
	b->k = makecomp(tfinalturn, tfinalang);	/* final turn ratio k */
	return b;
}

void branchcalc(struct branch *b) {
	struct comp adm; 
	struct comp k1;	/* k one */

	switch (b->type) {
	case AC: case FT:	/* ac line 暂时把FT也归入此类 */
		b->adm_se = b->adm_line;
		b->iadm_sh = makecomp(.0, b->linechar/2.0);
		b->jadm_sh = makecomp(.0, b->linechar/2.0);
		break;
	case 99:	/* fixed voltage ratio or fiexd phase angle */
		k1 = makecomp(1.0, .0);
		k1 = compmns(b->k, k1);			/* k - 1 */
		b->adm_se = compmul(b->adm_line, b->k); 		/* yk */
		adm = makecomp(.0, b->linechar/2.0);		/* 1/2 line char y0 */
		adm = compmul(adm, compmul(b->k ,b->k));	/* y0 . k^2 */
		adm = compadd(adm, 				/* + y . k . (k-1) */
			compmul(b->adm_line, compmul(b->k, k1)));
		b->iadm_sh = adm;
		adm = makecomp(.0, linechar/2.0);		/* y0 */
		adm = compadd(adm, compmul(b->adm_line, compinv(k1)));	/* +(1-k)y */
		b->jadm_sh = adm;
		break;
	default:
		fprintf(stderr, "cdf error: wrong branch type, %d\n", b->type);
		break;
	}
}

void nodecalc(struct node *t) {
	struct nodechain *pc;
	struct comp adm;

	adm = t->adm_sh;
	t->nconnect = 0;
	for (pc = t->nbr; pc != NULL; pc = pc->next) {
		adm = compadd(pc->b->adm_se, adm);
		if (t == pc->b->inode)
			adm = compadd(pc->b->iadm_sh, adm);
		else
			adm = compadd(pc->b->jadm_sh, adm);
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
