#include "ss.c"
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

/* determine if vsnprintf is available */
#ifndef HAVE_VSNPRINTF
#if ( defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__xlC__) || defined(_MSC_VER) )
#define HAVE_VSNPRINTF 1
#ifdef _MSC_VER
#ifndef vsnprintf
#define vsnprintf _vsnprintf
#endif /* VC define vsnprintf */
#endif /* _MSC_VER */
#endif /* good compilers */
#endif /* HAS_VSNPRINTF */

/* handle a trace command inside a format string 
 * return nonzero to quit and enter the normal mode */
static int tracecmd_(const char *cmd, va_list args,
    char **pfname, int *freq, int *verbose)
{
  const char *p;

  if ((p = strchr(cmd, '=')) == NULL)
    return 1;
  if (strncmp(cmd, "filename", 8) == 0) {
    p = va_arg(args, const char *);
    if (p != NULL) {
      sscpy(*pfname, p);
      if (*verbose)
        fprintf(stderr, "The trace file is now %s\n", *pfname);
    }
  } else if (strncmp(cmd, "freq", 4) == 0) {
    *freq = va_arg(args, int);
    if (*verbose)
      fprintf(stderr, "change frequency of writing trace to %d\n", *freq);
  } else if (strncmp(cmd, "verbose", 7) == 0) {
    *verbose = va_arg(args, int);
  } else {
    fprintf(stderr, "unknown command: %s\n", cmd);
    return 1;
  }
  return 0;
}

/* write trace in memory-buffered mode */
static int wtrace_low_(int cnt, int freq,
    const char *fname, const char *fmt, va_list args)
{
  const  int   maxmsg = 1024;
  static char *msg = NULL, *buf = NULL;
  int i;
  
  /* buffered mode */
  if (cnt == 0) { /* allocate msg and buf first */
    msg = ssnew(maxmsg * 2); /* leave some margin */
    buf = ssnew(maxmsg * freq);
  }
  if (msg == NULL || buf == NULL) {
    fprintf(stderr, "no buffer found.\n");
    return 1;
  }
  if (fmt != NULL) {
    if (strlen(fmt) >= (size_t) maxmsg) {
      fprintf(stderr, "the format string is too long.\n");
      return 1;
    }
#ifdef HAVE_VSNPRINTF 
    i = vsnprintf(msg, maxmsg, fmt, args);
#else /* we use vsprintf if vsnprintf is unavailable */
    i = vsprintf(msg, fmt, args);
#endif
    if (i >= maxmsg) {
      fprintf(stderr, "the message is too long.\n");
      return 1;
    }
    sscat(buf, msg);
  }

  /* flush buffered content to file, and possibly finish up */
  if ((cnt + 1) % freq == 0 || fmt == NULL) {
    FILE *fp;
    static const char *mode = NULL;

    if (buf && buf[0] != '\0') { /* in case nothing was written */
      xfopen(fp, fname, mode = (mode ? "a" : "w"), return -1);
      fputs(buf, fp);
      buf[0] = '\0';
      fclose(fp);
    }
    if (fmt == NULL) { /* finishing up */
      if (msg != NULL) ssdelete(msg); 
      if (buf != NULL) ssdelete(buf);
    }
  }
  return 0;
}

/* 
 * write trace using 
 * unbuffered (flags = 0) or buffered (flags = 1) output
 * but mixing the two modes is a major mistake!
 *
 * pass NULL to `fmt' to finish up
 * In addition, we support simple commands in `fmt'
 * If `fmt' starts with %@, we start the command mode
 * if fmt="%@filename=%s", the next argument is "tr.dat",
 *    then the trace file becomes tr.dat
 * if fmt="%@freq=%d", the next argument is 100,
 *    then we refresh every 100 calls.
 * */
int wtrace(const char *fmt, ...)
{
  static int   verbose = 1;
  static int   freq = 1000;
  static int   cnt = 0;
  static char *fname = NULL;
  va_list args;
  int i;

  if (fname == NULL) /* set the default file name */ 
    fname = ssdup("TRACE");

  /* start the command mode if the format string start with "%@"
   * the command mode allows setting parameters */
  if (fmt != NULL && fmt[0] == '%' && fmt[1] == '@') {
    if (cnt > 0) {
      fprintf(stderr, "trace: changing setting after something\n");
      return -1;
    }
    va_start(args, fmt);
    i = tracecmd_(fmt + 2, args, &fname, &freq, &verbose);
    va_end(args);
    if (i == 0) return 0;
  }

  va_start(args, fmt);
  i = wtrace_low_(cnt, freq, fname, fmt, args);
  va_end(args);

  /* fmt == NULL means finishing up, once = 0 for a fresh start */
  if (fmt == NULL) { /* finishing up */
    cnt  = 0;
    if (fname) ssdelete(fname);
  } else {
    cnt++;
  }
  return i;
}



#endif

