/* simplex.c */

#define N
#define M

Elm cm[] = {
	2, 3, 1, 0, 0, 0,
	3, 2, 0, 1, 0, 0,
	6, 1, 0, 0, 1, 0,
	1, 6, 0, 0, 0, 1 };
Elm d[] = {
	0, 0, 12, 12, 18, 18 };
Elm cost[] = {
	1, 2, 0, 0, 0, 0 };
Size idxcol[] = {
	0, 1, 2, 3, 4, 5 };

simplex(int lim, Elm eps) {
	Size i, j;

	
}

Size pos(Size i, Size j, Size ncol) {
	return i * ncol + j;
}

Elm getbase(Elm *v, Size i, Size j) {
	j += N-M; 
	return v[pos(i, idxcol[j], N)];
}

Elm getnbase(Elm *v, Size i, Size j) {
	return v[pos(i, idxcol[j], N)];
}

void setbase(Elm val, Elm *v, Size i, Sise j) {
	j += N-M;
	v[pos(i, idxcol[j], N)] = val;
} 
