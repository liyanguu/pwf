/* wrapper.h - KLU 皮函数头文件
 * 2017-5-17 */

#ifndef HEADER_wrapper
#define HEADER_wrapper

#include <cs.h>
#include <klu.h>

void wr_klu_defaults(void) ;
void wr_klu_analyze(int n, int *ap, int *ai) ;
void wr_klu_factor(int *ap, int *ai, double *ax) ;
int wr_klu_solve(int nldim, int nrhs, double *sol) ;

#endif
