#ifndef DIE_C__
#define DIE_C__

#include "die.h"

/* print a fatal message */
static void va_perr_(const char *file, int line, const char *why,
    const char *fmt, va_list args)
{
  fprintf(stderr, "error: ");
  vfprintf(stderr, fmt, args);
  if (fmt[strlen(fmt) - 1] != '\n')
    fprintf(stderr, "\n"); /* add a new line if needed */

  if (file != NULL) 
    fprintf(stderr, "file:  %s\n", file);
  if (line > 0)
    fprintf(stderr, "line:  %d\n", line);
  if (why != NULL) 
    fprintf(stderr, "why:   %s\n", why);
}

#define va_pmsg_(f, l, why, fmt, args)  vfprintf(stderr, fmt, args)

#define PERRMSG__(f, l, why, p, die)  \
    va_list args;                     \
    if (cond) {                       \
      va_start(args, fmt);            \
      p(f, l, why, fmt, args);        \
      va_end(args);                   \
      if (die) exit(1);               \
    }

/* print an message and die if cond is true */
#ifdef ZCHAVEVAM
void die_if_(const char *fsrc, int lnum, const char *why, 
    int die, int cond, const char *fmt, ...) 
{
  if (die) {
    PERRMSG__(fsrc, lnum, why, va_perr_, 1);
  } else {
    PERRMSG__(fsrc, lnum, why, va_pmsg_, 0);
  }
}
#else
void die_if(int cond, const char *fmt, ...)
{
  PERRMSG__(NULL, 0, NULL, va_perr_, 1);
}
void msg_if(int cond, const char *fmt, ...)
{
  PERRMSG__(NULL, 0, NULL, va_pmsg_, 0);
}
void fatal(const char *fmt, ...)
{
  int cond = 1;
  PERRMSG__(NULL, 0, NULL, va_perr_, 1);
}
#endif

#undef va_pmsg_
#undef PERRMSG__

#endif 

