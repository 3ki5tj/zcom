#ifndef DIE_C__
#define DIE_C__

#include "die.h"

#ifndef ZCHAVEVAM
#define err_fsrc_  NULL
#define err_lnum_  0
#define err_why_   NULL
#endif

/* print an error message and die */
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

/* die if cond is true */
#ifdef ZCHAVEVAM
void die_if_(const char *err_fsrc_, int err_lnum_, const char *err_why_, 
    int cond, const char *fmt, ...)
#else
void die_if(int cond, const char *fmt, ...)
#endif
{
  va_list args;

  if (cond) {
    va_start(args, fmt);
    va_perr_(err_fsrc_, err_lnum_, err_why_, fmt, args);
    va_end(args);
    exit(1);
  }
}

#ifndef ZCHAVEVAM
#ifdef DIE_LEGACY
void fatal(const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  va_perr_(err_fsrc_, err_lnum_, NULL, fmt, args);
  va_end(args);
  exit(1);
}
#endif
#undef  err_fsrc_
#undef  err_lnum_
#endif


#endif

