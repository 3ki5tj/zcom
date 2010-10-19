#ifndef DIE_C__
#define DIE_C__

#include "die.h"

/* print an error message */
void perrmsg__(const char *file, int line, const char *why,
    int err, const char *fmt, va_list args)
{
  if (err) fprintf(stderr, "error: ");
  vfprintf(stderr, fmt, args);
  if (fmt[strlen(fmt) - 1] != '\n')
    fprintf(stderr, "\n"); /* add a new line if needed */
  if (err) {
    if (file != NULL) fprintf(stderr, "file: %s\n", file);
    if (line > 0) fprintf(stderr, "line: %d\n", line);
    if (why != NULL && strcmp(why, "1") != 0) 
      fprintf(stderr, "cond: %s\n", why);
  }
}

#endif 

