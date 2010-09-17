#define ZCHAVEVAM 1

#ifndef DIE_H__
#define DIE_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#ifdef ZCHAVEVAM
#define fatal(fmt, ...)  die_if(1, fmt, ## __VA_ARGS__)
#define die_if(cond, fmt, ...)  \
  die_if_(__FILE__, __LINE__, #cond, 1, cond, fmt, ## __VA_ARGS__)
#define msg_if(cond, fmt, ...)  \
  die_if_(__FILE__, __LINE__, #cond, 0, cond, fmt, ## __VA_ARGS__)
void die_if_(const char *file, int line, const char *why, 
             int die, int cond, const char *fmt, ...);
#else
void die_if(int cond, const char *fmt, ...);
void msg_if(int cond, const char *fmt, ...);
void fatal(const char *fmt, ...);
#endif /* ZCHAVEVAM */

#endif

