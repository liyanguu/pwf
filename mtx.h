#ifndef H_MTX
#define H_MTX

#define ZEROEPS 1.0E-20
#define swapt(t, a, b) {t _z = a; a = b; b = _z;}
#define min(a, b) ((a) < (b) ? (a) : (b))
#define trymtxalloc(m, nrow, ncol) ((m==NULL) ? mtxalloc(nrow, ncol) : m)
typedef enum mtype_t { ROW_MAJOR, COL_MAJOR } mtype;

typedef int Size;
typedef double Elm; 
typedef Elm** Rv;
typedef struct mtx {
	Rv val;		/* ptrs to matrix elements' value */
	Size nrow;	/* # rows */
	Size ncol;	/* # cols */
	Size *ipiv;	/* index of pivot */
	mtype type;	/* matrix storage order */
} *Mtx;

Elm getel(Mtx, Size i, Size j);
void setel(Mtx, Size i, Size j, Elm);
void arraycpy(Elm *dst, Elm *src, Size lim);
void mtxread(Mtx, Elm *buf);
void readcol(Mtx, Size colno, Elm *buf);
void mtxwrite(Mtx, Elm *buf);
void mtxprint(char *title, Mtx m);
void xmtxprint(char *title, Mtx m);
void mtxnprt(char *title, ...);
Mtx mtxalloc(Size nrow, Size ncol);
void mtxfree(Mtx m);
Mtx mtxcopy(Mtx dst, Mtx src);
Mtx mtxdup(Mtx);
Mtx trans(Mtx);
Mtx exrow(Mtx m, Size row1, Size row2);
Mtx excol(Mtx m, Size c1, Size c2);
Mtx mtxmul(Mtx m, Mtx m1, Mtx m2);
Mtx xmtxmul(Mtx m, Mtx m1, Mtx m2);
Mtx mtxdif(Mtx m, Mtx m1, Mtx m2);
Mtx xmtxdif(Mtx m, Mtx m1, Mtx m2);
Mtx rowsub(Mtx m, Size *idx, Size nrow);
Mtx rowcopy(Mtx mc, Mtx m, Size *idx, Size nrow);
Mtx colsub(Mtx m, Size *idx, Size ncol);
Mtx colcopy(Mtx mc, Mtx m, Size *idx, Size ncol);
Mtx reshape(Mtx m, Size newrow, Size newcol);
int mtxsv(Mtx main, Mtx rhs);
Mtx mtxinv(Mtx m);
int ludcmp(Mtx, Size *, Elm *det);
void lubksb(Mtx, Size *, Elm *b);
void forward(Mtx a, Elm *b, Elm *y);
void back(Mtx, Elm *, Elm *);
void sparse(Elm *b, Size n, void (*asub)(Elm *, Elm *),
	void (*atsub)(Elm*, Elm*), Elm *x, Elm *rsq);
Mtx inverse(Mtx);
Size norm(Elm *v, Size lim, Elm *nr);

#endif
