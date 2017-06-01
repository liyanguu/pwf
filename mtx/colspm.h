/* 文件名：colspm.h (COLumn compressed SParse Matrix)
 * 用途： 稀疏矩阵的结构体定义
 * 日志：2017-5-23	修改与完善
	 2017-5-30      改为压缩列式结构  */

#ifndef HEADER_colspm
#define HEADER_colspm

#include "mtx.h"	/* 基本的量纲定义 */
#include "msg.h"

#define NONZ_BUFSZ 512

typedef struct colspm {
    	Elm *elms;	/* 矩阵的元素值 */
	Size *pos;	/* 行下标的位置索引 */
	Size *rows;	/* 行下标索引 */
	Size nz;	/* 已分配的非零元缓存大小 */
	Size dim;	/* 矩阵维数, 仅处理方阵 */
	Elm *colbuf;	/* 列元素缓存 */
	Size *index;	/* 行下标缓存 */
	Size colcnt;	/* 列元素计数 */
	Size curcol;	/* 当前的列，工作列 */
	Size cnt;	/* 总(非零)元素计数 */
} ColSpm;

ColSpm *colspm_alloc(void) ;
ColSpm *colspm_init(Size dim);
ColSpm *colspm_del(ColSpm *matx) ;
ColSpm *colspm_clr(ColSpm *matx) ;
int colspm_addcol(ColSpm *matx, Size colno) ;
int colspm_add(ColSpm *matx, Size rowno, Elm elm);
int colspm_endcol(ColSpm *matx);
void colspm_adjust(ColSpm *matx);
void colspm_print(ColSpm *);
void colspm_printex(ColSpm *matx);

#endif
