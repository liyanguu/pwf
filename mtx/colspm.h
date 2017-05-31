/* 文件名：colspm.h (COLumn compressed SParse Matrix)
 * 用途： 稀疏矩阵的结构体定义
 * 日志：2017-5-23	修改与完善
	 2017-5-30      改为压缩列式结构  */

#ifndef HEADER_colspm
#define HEADER_colspm

#include "mtx.h"	/* 基本的量纲定义 */
#include "msg.h"

typedef struct colspm {
    	Elm *elms;	/* 矩阵的元素值 */
	Size *pos;	/* 行下标的位置索引 */
	Size *rows;	/* 行下标索引 */
	Size nonz;	/* 非零元数 */
	Size dim;	/* 矩阵维数, 仅处理方阵 */
	Elm *colbuf;	/* 列元素缓存 */
	Size colcnt;	/* 列元素计数 */
	Size curcol;	/* 当前的列，工作列 */
	Size cnt;	/* 总元素计数 */
} ColSpm;

ColSpm *colspm_alloc(void) ;
ColSpm *colspm_init(Size dim, Size nonz) ;
ColSpm *colspm_del(ColSpm *matx) ;
ColSpm *colspm_clr(ColSpm *matx) ;
int colspm_addcol(ColSpm *matx, Size colno) ;
int colspm_add(ColSpm *matx, Size rowno, Elm elm);
void colspm_endcol(ColSpm *matx);
void colspm_print(ColSpm *);
void colspm_printex(ColSpm *matx);

#endif
