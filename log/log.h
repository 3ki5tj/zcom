#include "util.h"
#ifndef LOG_H__
#define LOG_H__



typedef struct {
  FILE *fp;
  const char *fname;
  unsigned flags;
} logfile_t;

#define LOG_WRITESCREEN  0x01
#define LOG_FLUSHAFTER   0x02
#define LOG_NOWRITEFILE  0x10
#define LOG_APPEND       0x80



INLINE logfile_t *log_open(const char *fn)
{
  logfile_t *log;

  xnew(log, 1);
  if (fn == NULL) fn = "LOG";
  log->fname = fn;
  log->flags = 0;
  return log;
}



INLINE int log_printf(logfile_t *log, char *fmt, ...)
{
  va_list args;

  if (log == NULL) return 1;
  if (log->fp == NULL) {
    const char *aw = (log->flags & LOG_APPEND) ? "a" : "w";
    xfopen(log->fp, log->fname, aw, return 1);
  }
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



INLINE void log_close(logfile_t *log)
{
  if (log == NULL) return;
  if (log->fp != NULL) { fclose(log->fp); log->fp = NULL; }
  free(log);
}



/* close & reopen log file to make sure that stuff is written to disk */
INLINE int log_hardflush(logfile_t *log)
{
  if (log->fp == NULL || log->fname == NULL) return 1;
  fclose(log->fp);
  xfopen(log->fp, log->fname, "a", return 1);
  return 0;
}

#if defined(HAVEVAM) && defined(NEED_WTRACE)
logfile_t log_stock_[1] = {{ NULL, "TRACE", 0 }};
#define wtrace(fmt, ...) { \
  if (fmt) log_printf(log_stock_, fmt, ##__VA_ARGS__); \
  else if (log_stock_->fp) { fclose(log_stock_->fp); log_stock_->fname = NULL; } }
#endif

#endif

