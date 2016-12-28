#ifndef MSG_H
#define MSG_H

#include <stdio.h>

extern int print_test_msg;

void msg(FILE *, char *, ...);
void exclude(char *);

#endif
