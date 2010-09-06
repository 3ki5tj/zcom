#ifndef CFG_H__
#define CFG_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct {
  int    n;           /* number of lines */
  int    canfree;     /* whether the struct is dynamically allocated */
  char **key,**value; /* key[i] = value[i] */
  char  *buf;         /* the whole configuration file */
} cfgdata_t;

cfgdata_t *cfgopen(const char *filenm);
void cfgclose(cfgdata_t *cfg);
int cfgget(cfgdata_t *cfg, void *var, const char *key, const char* fmt);

/* 
 * Helper macros for automating reading variables
 */

/* return a format string from a type
 * undetermined format goto "%%"  */
#define CFG_TP2FMT_(tp) \
    (strcmp(#tp, "int")       == 0 ? "%d"  : \
     strcmp(#tp, "unsigned")  == 0 ? "%u"  : \
     strcmp(#tp, "float")     == 0 ? "%f"  : \
     strcmp(#tp, "double")    == 0 ? "%lf" : \
     strcmp(#tp, "char *")    == 0 ? "%s"  : \
     strcmp(#tp, "null")      == 0 ? ""    : "%%")

#define CFG_PRINT_FILE_LINE_() \
    fprintf(stderr, "file: %s, line: %d\n", __FILE__, __LINE__); 

#ifndef CFG_FATAL_ACTION_  /* for programmer's mistake */
#define CFG_FATAL_ACTION_  exit(1)
#endif

/* print a message leading by `msg', and execute `action' 
 * low level error handling */
#define CFG_ERRMSG_(var, key, tp, def, desc, msg, action) {           \
    fprintf(stderr, "%s: var: %s, key: %s, type: %s\n",               \
        msg, #var, key, #tp);                                         \
    if (desc != NULL) fprintf(stderr, "desc: %s, ", desc);            \
    CFG_PRINT_FILE_LINE_()                                            \
    action; }

/* first guess a format, then read `var' through `name',
 * if the variable is missing, mismsg is printed, with misact
 * is called if the entry is not found in the configuration file   
 * the whole process is skipped if null is passed to tp 
 * */
#define CFG_GET_(cfg, var, name, tp, def, desc, mismsg, misact) {     \
  char *fmt_ = CFG_TP2FMT_(tp);                                       \
  if (fmt_[0] == '\0') {                                              \
    /* do nothing */ ;                                                \
  } else if (fmt_[1] == '%') {  /* unable to determine format */      \
    fprintf(stderr, "cannot determine format for %s\n", #tp);         \
    CFG_PRINT_FILE_LINE_()                                            \
    CFG_FATAL_ACTION_; /* fatal: programmer's mistake */              \
  } else if (sizeof(var) != sizeof(tp)) {                             \
    fprintf(stderr, "var. %s is not of type %s\n", #var, #tp);        \
    CFG_PRINT_FILE_LINE_()                                            \
    CFG_FATAL_ACTION_; /* fatal: programmer's mistake */              \
  } else if (0 != cfgget(cfg, &(var), name, fmt_)) {                  \
    CFG_ERRMSG_(var, name, tp, def, desc, mismsg, misact);            \
  }  }

/* conditionally (t0) get var from configuration file, 
 * then check if condition t1 is met, if not, misact is executed
 * neither t0 nor t1 should be empty, but can be 1 or 0;
 * finally, expression `eval', which can be empty, is evaluated;
 * a default value 'def' is always assigned at the beginning for safety,
 * the design is based on the concern that if anything fails
 * an uninitialized variable can be dangerous
 * */
#define CFG_GETC_(cfg, var, key, tp, def, t0, t1, eval, desc, mismsg, misact) \
  var = def;                                                                  \
  if (CFG_TEST_(t0)) {                                                        \
    CFG_GET_(cfg, var, key, tp, def, desc, mismsg, misact)                    \
    CFG_TESTERR_(t1, misact)                                                  \
    eval;                                                                     \
  }

/* test if expression t is true
 * t should *not* be empty 
 * t can be 1; since we tested #t[0], compilers usually don't complain
 *   the condition is always true; a further safer choice is probably
 *   strcmp(#t, "1") == 0
 * */
#define CFG_TEST_(t) ((#t[0] == '1' && #t[1] == '\0') ? 1 : (t))

/* if condition 'cond' is not true, execute `erract' */
#define CFG_TESTERR_(t, erract)                                       \
  if ( !CFG_TEST_(t) )  {                                             \
    fprintf(stderr, "failed cond: %s\n", #t);                         \
    CFG_PRINT_FILE_LINE_()                                            \
    erract;                                                           \
  }


/* allocate array, and initialize it with function def()
 * tp should be a pointer type, like int *, or double *,
 * def can be a number, or an expression using the index ia_  */
#define CFG_DARR_(arr, tp, def, cnt, desc, errmsg, erract)                  \
  if (sizeof(arr[0]) != sizeof(*((tp)arr))) {                               \
    fprintf(stderr, "array %s is not of type %s\n", #arr, #tp);             \
    CFG_ERRMSG_(arr, #arr, tp, def, desc, "type error", CFG_FATAL_ACTION_); \
  } else if ((arr = calloc( (cnt), sizeof(arr[0]) )) == NULL) {             \
    CFG_ERRMSG_(arr, #arr, tp, def, desc, errmsg, erract);                  \
  } else {                                                                  \
    size_t ia_;                                                             \
    for (ia_ = 0; ia_ < (cnt); ia_++) arr[ia_] = (def);                     \
  }


#endif  /* CFG_H__ */




/* print a mismatch message and die */
#define CFG_MISMATCH_TPFMT_(fmt, tp) \
    fprintf(stderr, "format %s doesn\'t match type %s", fmt, #tp);    \
    exit(1);

/* check if a type matches format by sizeof */
#define CFG_CHECKTPFMT(tp, fmt)                                       \
  if (fmt[0] != '%') {                                                \
    fprintf(stderr, "invalid format string %s\n", fmt);               \
    exit(1);                                                          \
  } else if (strchr("di", fmt[1]) && strcmp(#tp, "int") != 0) {       \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  } else if (strchr("uxo", fmt[1]) && strcmp(#tp, "unsigned") != 0) { \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  } else if (strchr("efg", fmt[1]) && strcmp(#tp, "float") != 0) {    \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  } else if (strcmp(fmt, "%lf") == 0 && strcmp(#tp, "double") != 0) { \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  } else if (strchr("sc", fmt[1]) && strcmp(#tp, "char *") != 0) {    \
    CFG_MISMATCH_TPFMT_(fmt, tp);                                     \
  }

/* read variable var in format fmt */
#define CFG_GETATOM0_(var, name, tp, fmt) {       \
  CFG_CHECKTPFMT(tp, fmt)                         \
  err = cfgget(cfg, &(var), name, fmt); }


