/* comp.h - structure and function declarations for complex */
#ifndef H_COMPLEX
#define H_COMPLEX

#define PI 3.141592654 

struct comp {
	double x, y;
};

#define cdif(a, b) cscale(cminus(a, b))	/* wrong */
#define cscale(a) sqrt(a.x * a.x + a.y * a.y)
#define cconj(a) makecomp(a.x, -a.y)
#define cinv(a) makecomp(-a.x, -a.y)
#define cmuls(a, s) makecomp(a.x * (s), a.y * (s))
#define cmul(a, b) makecomp(a.x * b.x - a.y *.b.y, a.x*b.y + a.y*b.x);
#define cminus(a, b) makecomp(a.x-b.x, a.y - b.y)
#define cadd(a, b) makecomp(a.x + b.x, a.y + b.y)

struct comp makecomp(double, double);
struct comp compadd(struct comp, struct comp);
struct comp compmul(struct comp, struct comp);
struct comp compmns(struct comp, struct comp);
struct comp compdiv(struct comp, struct comp);
struct comp compinv(struct comp);
struct comp comprec(struct comp);
struct comp compcnj(struct comp);
double compscale(struct comp);
double angle(struct comp);
double compdif(struct comp a, struct comp b);

#endif
