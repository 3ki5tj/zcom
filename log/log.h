#ifndef LOG_H__
#define LOG_H__

#include <stdio.h>
#include <stdarg.h>

typedef struct {
  FILE *fp;
  char *fname;
  int flag;
} logfile_t;

#define LOG_WRITESCREEN  0x01
#define LOG_FLUSHAFTER   0x02
#define LOG_NOWRITEFILE  0x10

logfile_t *logopen(char *filenm);
int logprintf(logfile_t *log, char *fmt, ...);
int loghardflush(logfile_t *log);
void logclose(logfile_t *log);

#endif

