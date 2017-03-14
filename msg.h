#ifndef MSG_H
#define MSG_H

#include <stdio.h>

extern int print_test_msg;

void msg(FILE *, char *, ...);
void excluding(char *);
void including(char *name);
int inlist(char *name, char **list);

#endif
