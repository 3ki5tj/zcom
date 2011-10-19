#ifndef INLINE
#define INLINE __inline static
#endif
#ifndef LOG_H__
#define LOG_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

typedef struct {
  FILE *fp;
  const char *fname;
  unsigned flags;
} logfile_t;

#define LOG_WRITESCREEN  0x01
#define LOG_FLUSHAFTER   0x02
#define LOG_NOWRITEFILE  0x10

logfile_t *log_open(const char *filenm);
int log_printf(logfile_t *log, char *fmt, ...);
void log_close(logfile_t *log);

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

