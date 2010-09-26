#ifndef DIE_H__
#define DIE_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#ifdef DIE_LEGACY__
void die_if_(const char *f, int l, int cond, const char *fmt, ...);
void msg_if_(const char *f, int l, int cond, const char *fmt, ...);
void fatal_(const char *f, int l, const char *fmt, ...);
#endif

void die_if(int cond, const char *fmt, ...);
#ifdef USE_MSG_IF
void msg_if(int cond, const char *fmt, ...);
#endif
#ifdef USE_FATAL
void fatal(const char *fmt, ...);
#endif

#endif

