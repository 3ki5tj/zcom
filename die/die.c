#ifndef DIE_C__
#define DIE_C__

#include "die.h"

/* print a fatal message */
static void va_perr_(const char *file, int line, int err, 
    const char *fmt, va_list args)
{
  if (err) fprintf(stderr, "error: ");
  vfprintf(stderr, fmt, args);
  if (fmt[strlen(fmt) - 1] != '\n')
    fprintf(stderr, "\n"); /* add a new line if needed */

  if (err) {
    if (file != NULL) fprintf(stderr, "file:  %s\n", file);
    if (line > 0) fprintf(stderr, "line:  %d\n", line);
  }
}

#define PERRMSG__(f, l, c, x) {   \
  va_list args;                   \
  if (#c[0] == '1' || c) {        \
    va_start(args, fmt);          \
    va_perr_(f, l, x, fmt, args); \
    va_end(args);                 \
    if (#x[0] == '1') exit(1);    \
  } }

#ifdef DIE_LEGACY__
/* clumsy version */
void die_if_(const char *f, int l, int cond, const char *fmt, ...) 
  PERRMSG__(f, l, cond, 1)
void msg_if_(const char *f, int l, int cond, const char *fmt, ...) 
  PERRMSG__(f, l, cond, 0)
void fatal_(const char *f, int l, const char *fmt, ...) 
  PERRMSG__(f, l, 1, 1)
#endif

/* shorter version */
void die_if(int cond, const char *fmt, ...) PERRMSG__(NULL, 0, cond, 1)
#ifdef USE_MSG_IF
void msg_if(int cond, const char *fmt, ...) PERRMSG__(NULL, 0, cond, 0)
#endif
#ifdef USE_FATAL
void fatal(const char *fmt, ...) PERRMSG__(NULL, 0, 1, 1)
#endif

#undef PERRMSG__

#endif 

