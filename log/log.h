#ifndef LOG_H__
#define LOG_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

typedef struct {
  FILE *fp;
  const char *fname;
  unsigned flag;
} logfile_t;

#define LOG_WRITESCREEN  0x01
#define LOG_FLUSHAFTER   0x02
#define LOG_NOWRITEFILE  0x10

logfile_t *log_open(const char *filenm);
int log_printf(logfile_t *log, char *fmt, ...);
int log_hardflush(logfile_t *log);
void log_close(logfile_t *log);

#ifdef HAVEVAM
logfile_t log_stock_[1] = {{ NULL, "TRACE", 0 }};
#define wtrace(fmt, ...) { \
  if (fmt) log_printf(log_stock_, fmt, ##__VA_ARGS__); \
  else if (log_stock_->fp) { fclose(log_stock_->fp); log_stock_->fname = NULL; } } 
#endif

#endif

