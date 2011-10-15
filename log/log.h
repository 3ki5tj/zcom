#ifndef LOG_H__
#define LOG_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

typedef struct {
  FILE *fp;
  const char *fname;
  int flag;
} logfile_t;

#define LOG_WRITESCREEN  0x01
#define LOG_FLUSHAFTER   0x02
#define LOG_NOWRITEFILE  0x10

logfile_t *log_open(const char *filenm);
int log_printf(logfile_t *log, char *fmt, ...);
int log_hardflush(logfile_t *log);
void log_close(logfile_t *log);

#endif

