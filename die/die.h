#define ZCINLINE __inline static

#ifndef DIE_H__
#define DIE_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

void perrmsg__(const char *file, int line, const char *why,
    int err, const char *fmt, va_list args);

#ifdef ZCHAVEVAM

ZCINLINE void perrmsg_(const char *file, int line, const char *why,
    int cond, int err, const char *fmt, ...)
{
  va_list args;
  
  if (cond) {
    va_start(args, fmt);
    perrmsg__(file, line, why, err, fmt, args);
    va_end(args);
    if (err) exit(1);
  }
}
#define die_if(cond, fmt, ...) \
  perrmsg_(__FILE__, __LINE__, #cond, cond, 1, fmt, ## __VA_ARGS__)
#define msg_if(cond, fmt, ...) \
  perrmsg_(__FILE__, __LINE__, #cond, cond, 0, fmt, ## __VA_ARGS__)
#define fatal(fmt, ...)  die_if(1, fmt, ## __VA_ARGS__)

#else /* ZCHAVEVAM */

#define PERRMSG__(c, x) {                     \
  va_list args;                               \
  if ((#c[0] == '1' && #c[1] == '\0') || c) { \
    va_start(args, fmt);                      \
    perrmsg__(NULL, -1, #c, x, fmt, args);    \
    va_end(args);                             \
    if (#x[0] == '1') exit(1);                \
  } }
void die_if(int cond, const char *fmt, ...) PERRMSG__(cond, 1)
#ifdef USE_MSG_IF
void msg_if(int cond, const char *fmt, ...) PERRMSG__(cond, 0)
#endif
#ifdef USE_FATAL
void fatal(const char *fmt, ...) PERRMSG__(1, 1)
#endif
#undef PERRMSG__

#endif /* ZCHAVEVAM */


#endif

