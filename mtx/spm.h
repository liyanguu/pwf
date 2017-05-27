/* 文件名：spm.h(SParse Matrix)
 * 用途： 稀疏矩阵的结构体定义
 * 日志：2017-5-23	修改与完善 */

#ifndef Header_Spm
#define Header_Spm

#include "mtx.h"	/* 基本的量纲定义 */

typedef struct spel_struct {		/* a sparse matrix element */
	Elm v;		/* value */
	Size colno;	/* coloumn no, start from 0 */
	Size rowno;	/* row no, start from 0 */
	struct spel_struct *colnext;	/* ptr to the next col element */
	struct spel_struct *rownext;	/* ptr to the next row element */
} *Spel; 

typedef struct spm_struct {		/* a sparse matirx */
	Spel *colhead;		/* ptr to cols */ 
	Spel *rowhead;		/* ptr to rows */ 
	Size ncol;	/* # coloumns */
	Size nrow;	/* # rows */
} *Spm;

typedef struct compressed_spm_struct {
	Size nnz;	/* # non-zero elements */
	Size *cnz;	/* # non-zero elements in a column */
	Size *cnz_cnt;	/* a counter for cnz */
	Size *ap;	/* compressed-column data */
	Size *ai;
	Elm *ax;
} *CSpm;

#define validate(s, i, j) ((i) >= 0 && (i) < s->nrow && (j) >=0 && (j) < s->ncol)

void spprint(char *title, Spm a);
void spdelet(Spm a);
void sprem(Spm a, Size rowno, Size colno);
int spadd(Spm a, Size rowno, Size colno, Elm v);
int spget(Spm a, Size rowno, Size colno, Elm *v);
Spel spfind(Spm a, Size rowno, Size colno);
Spm spinit(Size nrow, Size ncol);
Spm spalloc(void);
Spel spelalloc(void);
void spcompress(Spm s, CSpm cs);

#endif
