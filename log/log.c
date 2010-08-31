#include "ss.h"

#ifndef LOG_C__
#define LOG_C__

/*
 * =======================================================================
 *
 * LOG file routines
 *
 * ========================================================================
 */
#include "log.h"

logfile_t *logopen(char *filenm)
{
  logfile_t *log;
  
  if (filenm == NULL) /* assign a default name */
    filenm = "LOG";
  if ((log=calloc(1, sizeof(*log))) == NULL) {
    fprintf(stderr, "cannot allocate memory for log file %s\n", filenm);
    return NULL;
  }
  /* We merely copy the name of the file,
   * the file is not opened until the first logprintf call */
  log->fname = ssdup(filenm);
  log->flag=0;
  return log;
}

int logprintf(logfile_t *log, char *fmt, ...)
{
  va_list args;

  if (log == NULL) return 1;

  if (log->fp == NULL) 
    log->fp = fopen(log->fname, "w");
  if (log->fp == NULL) {
    fprintf(stderr, "log [%s] cannot be opened.\n", log->fname);
    return 1;
  }
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
int loghardflush(logfile_t *log)
{
  if (log->fp == NULL || log->fname == NULL) 
    return 1;
  fclose(log->fp);
  if ((log->fp = fopen(log->fname, "a")) == NULL) {
    fprintf(stderr, "cannot reopen the log file [%s].\n",
        log->fname);
    return 1;
  }
  return 0;
}

void logclose(logfile_t *log) 
{
  if (log == NULL) 
    return;
  if (log->fp != NULL) {
    fclose(log->fp);
    log->fp = NULL;
  }
  if (log->fname != NULL) {
    ssdelete(log->fname);
  }
  free(log);
}

#endif

