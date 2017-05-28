/* memmange.c - 为数据结构分配与释放内存 
 * 2017-5-28	*** */

#include "pwf.h"

struct nodechain *chainalloc(void) {
	return malloc(sizeof(struct nodechain));
}

void chainfree(struct nodechain *p) {
    	if (p != NULL)
	    	free(p);
}

void nodechainfree(struct nodechain* p) {
	struct nodechain *q;

	for ( ; p != NULL; p = q) {
		q = p->next;
		chainfree(p);
	}
}

/* addchain: add a member to chain pc, with node pn, branch pb */
struct nodechain *addchain(struct nodechain *pc, struct node *pn, struct branch *pb) {
	struct nodechain *p;

	if ((p = chainalloc()) != NULL) {
		p->n = pn;
		p->b = pb;
		p->pw_f = makecomp(.0, .0);
		p->next = pc;
	}
	return p;
}

/* invchain: inverse a nodechain in place */
struct nodechain *invchain(struct nodechain *h) {
	struct nodechain *prev, *p;

	if (h == NULL)
		return NULL;
	for (prev = h; (p = h->next) != NULL; prev = p) {
		h->next = p->next;
		p->next = prev;
	}
	return prev;
}

struct node *nodealloc(void) {
	return malloc(sizeof(struct node));
}

void nodefree(struct node *t) {
	if (t != NULL) {
		free(t->name);
		nodechainfree(t->nbr);
		free(t);
	}
}

struct branch *branchalloc(void) {
	return malloc(sizeof (struct branch));
}

void branchfree(struct branch *b) {
	if (b != NULL)
		free(b);
}
