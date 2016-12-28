#ifndef H_MTX
#define H_MTX

#define ZEROEPS 1.0E-20
#define swapt(t, a, b) {t _z = a; a = b; b = _z;}

typedef int Size;
typedef double Elm;
typedef Elm *Rv;

typedef struct mtx {
	Rv *val;
	Size nrow, ncol;
} *Mtx;

struct spel {
	Elm v;
	Size colno;
	struct spel *nextcol;
};

typedef struct {
	struct spel *val;
	Size nrow, ncol;
} *Spm;


void mtxread(Mtx, Elm *);
void readcol(Mtx, Size, Elm *);
void mtxwrite(Mtx, Elm *);
void mtxcopy(Mtx dst, Mtx src);
void mtxprint(char *title, Mtx m);
void mtxnprt(int n, char *title, ...);
Mtx mtxalloc(Size nrow, Size ncol);
void mtxfree(Mtx m);
Mtx trans(Mtx);
Mtx mtxmul(Mtx m1, Mtx m2);
Mtx mtxdif(Mtx m1, Mtx m2);
Mtx rowsub(Mtx m, Size nrow, Size *idx);
Mtx colsub(Mtx m, Size ncol, Size *idx);
void exrow(Mtx m, Size row1, Size row2);
void excol(Mtx m, Size c1, Size c2);
int ludcmp(Mtx, Size *, Elm *det);
void lubksb(Mtx, Size *, Elm *b);
void forward(Mtx a, Elm *b, Elm *y);
void back(Mtx, Elm *, Elm *);
void sparse(Elm *b, Size n, void (*asub)(Elm *, Elm *),
	void (*atsub)(Elm*, Elm*), Elm *x, Elm *rsq);
Mtx inverse(Mtx);
Spm makesp(Size, Size);
void deletsp(Spm);
void printsp(char *, Spm);
void add(Spm, Size, Size, Elm);
void gauss_jordan(Spm, Elm *);

#endif
