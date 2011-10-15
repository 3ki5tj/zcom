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
  log->flag = 0;
  return log;
}

int log_printf(logfile_t *log, char *fmt, ...)
{
  va_list args;

  if (log == NULL) return 1;
  if (log->fp == NULL) xfopen(log->fp, log->fname, "w", return 1);
  if ((log->flag & LOG_NOWRITEFILE) == 0) {
    va_start(args, fmt);
    vfprintf(log->fp, fmt, args);
    va_end(args);
  }
  if (log->flag & LOG_WRITESCREEN) {
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
  if (log->flag & LOG_FLUSHAFTER)
    fflush(log->fp);
  return 0;
}

/* close & reopen log file to make sure that stuff is written to disk */
int log_hardflush(logfile_t *log)
{
  if (log->fp == NULL || log->fname == NULL) return 1;
  fclose(log->fp);
  xfopen(log->fp, log->fname, "a", return 1);
  return 0;
}

void log_close(logfile_t *log)
{
  if (log == NULL) return;
  if (log->fp != NULL) fclose(log->fp);
  free(log);
}

#endif

