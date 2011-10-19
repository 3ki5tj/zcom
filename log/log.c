#include "util.h"
#ifndef LOG_C__
#define LOG_C__
#include "log.h"

logfile_t *log_open(const char *fn)
{
  logfile_t *log;

  xnew(log, 1);
  if (fn == NULL) fn = "LOG";
  log->fname = fn;
  log->flags = 0;
  return log;
}

int log_printf(logfile_t *log, char *fmt, ...)
{
  va_list args;

  if (log == NULL) return 1;
  if (log->fp == NULL) xfopen(log->fp, log->fname, "w", return 1);
  if ((log->flags & LOG_NOWRITEFILE) == 0) {
    va_start(args, fmt);
    vfprintf(log->fp, fmt, args);
    va_end(args);
  }
  if (log->flags & LOG_WRITESCREEN) {
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
  if (log->flags & LOG_FLUSHAFTER)
    fflush(log->fp);
  return 0;
}

void log_close(logfile_t *log)
{
  if (log == NULL) return;
  if (log->fp != NULL) fclose(log->fp);
  free(log);
}

#endif

