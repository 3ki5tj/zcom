#define ZCHAVEVAM 1

#ifndef DIE_H__
#define DIE_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#ifdef ZCHAVEVAM
#define fatal(fmt, ...)         die_if(1, fmt, ## __VA_ARGS__)
#define die_if(cond, fmt, ...)  die_if_(__FILE__, __LINE__, #cond, \
                                        cond, fmt, ## __VA_ARGS__)
void die_if_(const char *file, int line, const char *why, 
             int cond, const char *fmt, ...);
#else
#ifdef DIE_LEGACY
void fatal(const char *fmt, ...);
#endif
void die_if(int cond, const char *fmt, ...);
#endif /* ZCHAVEVAM */

#endif

