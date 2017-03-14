/* comp.c - complex algorithms */
#include "comp.h"
#include <math.h>

/* makecomp: return a complex number (x, y) */
struct comp makecomp(double x, double y) {
	struct comp temp; 

	temp.x = x;
	temp.y = y;
	return temp;
}


/* compadd: add complex a and complex b, return the sum */
struct comp compadd(struct comp a, struct comp b) {
	a.x += b.x;
	a.y += b.y;
	return a;
}


/* compmns: minus complex a by complex b, return the result */
struct comp compmns(struct comp a, struct comp b) {
	a.x -= b.x;
	a.y -= b.y;
	return a;
}

/* compdiv: divide complex a by complex b, 
	if b = (0, 0) return (0, 0), else return the result */
struct comp compdiv(struct comp a, struct comp b) {
	struct comp temp;
	double buttom;
	
	buttom = b.x * b.x + b.y * b.y;
	if (buttom == 0)
		return makecomp(.0, .0);
	else {
		temp.x = (a.x * b.x + a.y * b.y) / buttom;
		temp.y = (a.y * b.x - a.x * b.y) / buttom;
		return temp;
	}
}


/* compmul: multiply complex a and complex b, return the result */
struct comp compmul(struct comp a, struct comp b) {
	struct comp temp;

	temp.x = a.x * b.x - a.y * b.y;
	temp.y = a.x * b.y + a.y * b.x;
	return temp;
}


/* compmuls: multiply complex a and scaler s, return the result */
struct comp compmuls(struct comp a, double s) {
	a.x *= s;
	a.y *= s;
	return a;
}


/* compinv: return the inverse of complex a = (x, y), i.e. (-x, -y) */
struct comp compinv(struct comp a) {
	a.x = - a.x;
	a.y = - a.y;
	return a;
}


/* compcnj: return the conjunction of complex a = (x, y), i.e. (x, -y) */
struct comp compcnj(struct comp a) {
	a.y = - a.y;
	return a;
}

/* comprec: return the reciprocal of complex a */
struct comp comprec(struct comp a) {
	struct comp temp;
	
	temp = makecomp(1.0, .0);
	return compdiv(temp, a);
}	


/* compscale: return the scale of complex a */
double compscale(struct comp a) {
	return sqrt(a.x * a.x + a.y * a.y);
}

/* angle: return the angle of complex a in rad, 
	from -PI to PI */
double angle(struct comp a) {
	if (fabs(a.x) == 1e-20) {
		if (fabs(a.y) == 1e-20)
			return 0;
		else if (a.y > 1e-20)
			return PI / 2;
		else
			return -PI /2;
	} else if (a.x > 1e-20) {
		if (fabs(a.y) == 1e-20)
			return 0;
		return atan(a.y / a.x);
	} else {
		if (fabs(a.y) == 1e-20)
			return PI;
		return atan(a.y / a.x) + ((a.y > 1e-20) ? PI : -PI); 
	}
}

/* compdif: return the difference of complex a and complex b,
	i.e. the scale of a - b */ 
double compdif(struct comp a, struct comp b) {
	return compscale(compmns(a, b));
}
