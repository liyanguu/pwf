/* 文件名：cdf.c
 * 用法：include "pwf.h"
 * 简述：读取或写出CDF格式文件的函数
   更改记录：
	2017-3-19	modify bus section data, add ctlmin,max
	2017-5-13~14	修改部分函数, 修正错误：不能正确读入CDF文件
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "pwf.h"
#include "comp.h"

		/* 共 6 项 */
#define TITLEFMT	" %8c %20c %lf %d %1s %30c"
#define TITLEPFMT	" %-8s %-20s %6.1f %4d %1s %-30s\n"

		/* 共 18 项 */
#define BUSFMT "%d %12c %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d"

#define BUSPFMT "%4d %12s %2d %2d %2d %6.3f%7.2f%9.1f%10.1f%8.1f%8.1f %7.1f %6.3f%8.1f%8.1f%8.1f%8.1f %4d\n"

		/* 共 21 项 */
#define BRANCHFMT "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf"

#define BRANCHPFMT "%4d %4d%3d%3d %1d %1d%10.5f%10.5f%10.5f%6d%6d%6d%4d %d %7.4g %7.4g%7.4g%7.4g%7.5g %7.4g%7.4g\n"


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
static char line[MAXLINE];	/* CDF 行缓存 */
static int section; 		/* 数据区标志 */
static int lineno = 0;		/* 行号 */

/* readcdfline: 读取CDF文件中的下一行，并去除注释 */
char *readcdfline(FILE *fp) {
    	char *ps;
	int trim(char *);

	while ((ps = fgets(line, MAXLINE, fp)) != NULL) {
	    	lineno++;
		if (ps[0] != '#') {
		    	trim(ps);
		    	break;
		}
	}
	return ps;
}

/* readcdf: 读取CDF文件fp，完成节点与支路的创建 */
int readcdf(FILE *fp) {
	char *ps;

	if (!titlescan(readcdfline(fp))) {
		fprintf(stderr, "readcdf: missing title card\n");
		return 0;
	}
	while ((ps = readcdfline(fp)) != NULL) {
		if(strstr(ps, "BUS DATA FOLLOWS") == ps) {
			section = BUS;
		} else if(strstr(ps, "BRANCH DATA FOLLOWS") == ps) {
			section = BRANCH;
		} else if (strstr(ps, "-999") == ps
				|| strstr(ps, "-99") == ps
				|| strstr(ps, "-9") == ps) {
			section = END;
		} else if (strstr(ps, "END OF DATA") == ps) {
			section = EOD;
		} else {
			switch (section) {
			case BUS:
				if (!busscan(ps)) {
					return 0;
				}
				break;
			case BRANCH:
				if (!branchscan(ps)) {
					return 0;
				}
				break;
			default:
				fprintf(stderr, "readcdf: wrong section %d at line %d\n", 
					section, lineno);
			}
		}
	}
	if (section != EOD) {
		fprintf(stderr, "readcdf: missing END OF DATA\n");
		return 0;
	}
	return 1;
}

/* writecdf: write load flow result to file fp in CDF FORMAT */
/* writecdf: 将潮流计算结果以CDF格式写入文件 fp 中 */
int writecdf(FILE *fp) {
	struct node *t;
	struct branch *b;
	int i, maxnode, maxbranch;

	getsize(&maxnode, &maxbranch);
	writetitle(line);
	fprintf(fp, "%s", line);

	fprintf(fp, "BUS DATA FOLLOWS %30d ITEMS\n", maxnode);
	for (i = 0; (t=getnode(i)) != NULL; i++) {
		writebus(line, t);
		fprintf(fp, "%s", line);
	}
	fprintf(fp, "-999\n");

	fprintf(fp, "BRANCH DATA FOLLOWS %30d ITEMS\n", maxbranch);
	for (i = 0; (b=getbranch(i)) != NULL; i++) {
		writebranch(line, b);
		fprintf(fp, "%s", line);
	}
	fprintf(fp, "-999\n");

	fprintf(fp, "END OF DATA\n");
	return !ferror(fp);
}

/* closecdf: 关闭文件fp，并将行号归零 */
void closecdf(FILE *fp) {
	fclose(fp);
	lineno = 0;
}

/* titlescan: 读入标题区的信息 */
int titlescan(char *linebuf) {
	int n;

	n = sscanf(linebuf, TITLEFMT, senddate, sendername, &basemva,
			&year, season, id);
	if (n != 6) {
		fprintf(stderr, 
		"titlescan: can't read title at line %d:\n%s\n", lineno, linebuf);
		return 0;
	} 
	senddate[8] = '\0';
	sendername[20] = '\0';
	season[1] = '\0';
	id[30] = '\0';
	return 1;
}

/* busscan: 读入母线信息 */
int busscan(char *linebuf) {
	struct node *pn;
	int n;

	if ((pn = addnode()) == NULL) {
		fprintf(stderr, "busscan: can't add new node\n");
		return 0;
	}
	n = sscanf(linebuf, BUSFMT, 
			&busno, name, 
			&lfano, &lzno, &bustype,
			&finalvolt, &finalang, 
			&loadmw, &loadmvar,
			&genmw, &genmvar, 
			&basekv, &volt_ctl, 
			&opmax, &opmin, 
			&shcon, &shsus, 
			&rcbusno);
	if (n != 18) {
		fprintf(stderr, 
		"busscan: can't read bus at line %d:\n%s\n", lineno, linebuf);
		return 0;
	}
	name[12] = '\0';
	makenode(pn);
	fprintf(stderr, "busscan: scaned bus %d\n", pn->no);
	return 1;
}

/* branchscan: 读入支路信息 */
int branchscan(char *linebuf) {
	struct branch *pb;
	int n;

	if ((pb = addbranch()) == NULL) {
		fprintf(stderr, "branchscan: can't input new branch\n");
		return 0;
	}
	n = sscanf(linebuf, BRANCHFMT, 
	&tbusno, &zbusno, 
	&lfano, &lzno,
	&circuit, &branchtype, 
	&br, &bx, &linechar,
	&mvarate1, &mvarate2, &mvarate3,
	&ctlbusno, &ctlside, 
	&tfinalturn, &tfinalang,
	&mintap, &maxtap, 
	&stepsize, 
	&minvolt, &maxvolt);
	if (n != 21) {
		fprintf(stderr, "branchscan: at cdf line %d:\n%s\n", lineno, linebuf);
		return 0;
	}
	fprintf(stderr, "branchscan: scaned branch %d-%d\n", tbusno, zbusno);
	return (makebranch(pb) == NULL) ? 0 : 1;
}

/* writetitle: 将标题信息写入行缓存linebuf中 */
void writetitle(char *linebuf) {
	time_t mytime;
	struct tm *loctime;

	time(&mytime);
	loctime = localtime(&mytime);
	strftime(senddate, sizeof(senddate), "%D", loctime);
	strncpy(sendername, "WHU LAB", sizeof(sendername)-1);
	sendername[20] = '\0';

	sprintf(linebuf, TITLEPFMT, senddate, sendername, basemva,
		year, season, id);
}

void writebus(char *linebuf, struct node *t) {
	double genmw, genmvar;
	double loadmw, loadmvar;
	double pmax, pmin;
	int bustype;
	
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
	loadmw = t->pload * basemva;
	loadmvar = t->qload * basemva;
	genmw = t->pw.x * basemva + loadmw;
	genmvar = t->pw.y * basemva + loadmvar;

	sprintf(linebuf, BUSPFMT,
		t->no, t->name, 1, 1, bustype,
		compscale(t->volt), angle(t->volt) * 180.0 / PI,
		fixzero(loadmw), fixzero(loadmvar),
		fixzero(genmw), fixzero(genmvar),
		t->basekv, t->volt_ctl,
		pmax, pmin,
		t->adm_sh.x, t->adm_sh.y, 0);
}

void writebranch(char *linebuf, struct branch *b) {
	struct comp imp;

	imp = comprec(b->adm_line);

	sprintf(linebuf, BRANCHPFMT,
		b->inode->no, b->jnode->no,
		1, 1, 1, b->type,
		imp.x, imp.y, b->linechar, 
		0, 0, 0, 0, 0, b->t, .0,
		.0, .0, .0, .0, .0);
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
/* makenode: 从已读入的CDF母线数据中创建一个节点 */
struct node *makenode(struct node *t) {
	double ang_rad;

	t->no = busno;
	t->flag = 0;
	t->type = bustype;
	t->name = strdup(name);
	t->basekv = basekv;
	t->nbr = NULL;
	t->pload = loadmw / basemva;
	t->qload = loadmvar / basemva;
	if (t->type == PV) {
	/* for PV nodes, opmin = min MVAr, opmax = max MVAr */
		t->q_min = opmin / basemva;
		t->q_max = opmax / basemva;
	} else if (t->type == PQ) {
	/* for PQ nodes, minop = min voltage limit, maxop = max volt limit */
		t->vt_min = opmin;
		t->vt_max = opmax;
		t->var_min = ctlmin;
		t->var_max = ctlmax;
	}
	t->volt_ctl = volt_ctl;
	t->adm_sh = makecomp(shcon, shsus);
	t->adm_self = t->adm_sh; /* 初始化节点的自导纳为其并联导纳 */ 
	ang_rad = finalang / 180.0 * 3.141592654;
	t->volt = makecomp(finalvolt * cos(ang_rad), 
			finalvolt * sin(ang_rad));
	t->pw = makecomp((genmw - loadmw)/basemva,
			(genmvar - loadmvar)/basemva);
	return t;
}

/* makebranch: 从已读入的支路信息创建一条支路, 
 * 建立节点-支路连接关系，计算其pi等效模型,
 * 并计算支路两端节点的自导纳 */
struct branch *makebranch(struct branch *b) {
	struct comp imp_b;
	struct node *ipn, *jpn;

	ipn = findnode(tbusno);	
	jpn = findnode(zbusno);
	if (ipn == NULL || jpn == NULL) {
		fprintf(stderr, 
		"makebranch: line without terminal nodes:\n%s\n", line);
		return NULL;
	}
	b->inode = ipn;
	b->jnode = jpn;
	b->type = branchtype;
	b->linechar = linechar;
	imp_b = makecomp(br, bx);
	b->adm_line = comprec(imp_b);
	b->t = tfinalturn;
	branchcalc(b); /* 计算pi等效模型 */
	/* 增加相邻节点 */
	ipn->nbr = addchain(ipn->nbr, b->jnode, b);	
	jpn->nbr = addchain(jpn->nbr, b->inode, b);
	/* 计算支路两端节点的自导纳 */
	nodebranch(ipn, b);
	nodebranch(jpn, b);
	return b;
}

/* calculate the pi-model of transmission lines & transformers 
   return the branch mutual admitance */
/* branchcalc: 计算线路或变压器的pi型模型,并返回支路的互导纳 */
struct comp branchcalc(struct branch *b) {
	struct comp adm; 
	double t;

	switch (b->type) {
	case AC:	/* ac line */
		b->adm_se = b->adm_line;
		b->adm_ish = b->adm_jsh = makecomp(.0, b->linechar/2.0);
		break;
	case FT: case VT:	/* fixed & variable tap transformers */
		adm = b->adm_line;
		t = b->t;
		b->adm_se = compmuls(adm, 1/t); /* Tij * yij */
		b->adm_ish = compmuls(adm, (1-t)/(t * t)); /* (Tij ** 2 - Tij)*yij*/
		b->adm_jsh = compmuls(adm, (t-1)/t);  /* (1 - Tij) * yij */
		break;
	default:
		fprintf(stderr, "cdf error: wrong branch type, %d\n", b->type);
		return makecomp(.0, .0); /* 返回复数零 */
	}
	return b->adm_mut = compinv(b->adm_se);
}

/* calculate and return the self-admitance of the node :
   Yii = yii0 + SUM OF (yik0 + yik) */
/* nodecalc: 计算并返回节点的自导纳:
 * 计算公式 Yii = yii0 + SUM OF (yik0 + yik) */
struct comp nodecalc(struct node *t) {
	struct nodechain *pc;
	struct comp adm;

	adm = t->adm_sh; /* yii0 */
	t->nconnect = 0;
	for (pc = t->nbr; pc != NULL; pc = pc->next) {
		adm = compadd(pc->b->adm_se, adm); /* yik */
		if (t == pc->b->inode)
			adm = compadd(pc->b->adm_ish, adm); /* yik0 */
		else
			adm = compadd(pc->b->adm_jsh, adm);
		t->nconnect++;
	}
	return t->adm_self = adm;
}

/* nodebranch: 由支路pb的pi型参数计算节点pn的自导纳 */
void nodebranch(struct node *pn, struct branch *pb) {
    	struct comp adm = pn->adm_self;

	/* yik0 */
   	if (pn == pb->inode) {
	    	adm = compadd(adm, pb->adm_ish);
	} else if (pn == pb->jnode) {
	    	adm = compadd(adm, pb->adm_jsh);
	} else {
	    	fprintf(stderr, "nodebrnch: node doesn't relate to branch\n");
		return;
	}
	adm = compadd(adm, pb->adm_se); /* yik */
	pn->adm_self = adm;
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
