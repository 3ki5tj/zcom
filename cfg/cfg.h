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
 * Helper macros for reading structure member
 * 
 * assume the following variables are defined
 * cfgdata_t *cfg;
 *
 * define macro CFG_STPTR as a pointer variable to the structure
 *
 * define macro CFG_STPFX as the common prefix for variables in the 
 * configuration file, e.g., if a configuration file has
 *   abc_color = yellow
 *   abc_number = 3
 * then CFG_STPFX is "abc_"
 *
 * a label ERR at the function exit
 *
 * */

/* return a format string from a type
 * undetermined format goto "%%"  */
#define CFG_TP2FMT(tp) \
    (strcmp(#tp, "int")       == 0 ? "%d"  : \
     strcmp(#tp, "unsigned")  == 0 ? "%u"  : \
     strcmp(#tp, "float")     == 0 ? "%f"  : \
     strcmp(#tp, "double")    == 0 ? "%lf" : \
     strcmp(#tp, "char *")    == 0 ? "%s"  : "%%")

#define CFG_MISERROR(var, key, tp, def)                               \
    fprintf(stderr, "cannot read variable %s (key: %s) of type %s, "  \
        "file: %s, line: %d\n", #var, key, #tp, __FILE__, __LINE__);  \
    goto ERR;

#define CFG_MISDFLT(var, key, tp, def)                                      \
    fprintf(stderr, "assume default value %s for %s (key: %s) of type %s, " \
        "file: %s, line: %d\n", #def, #var, key, #tp, __FILE__, __LINE__); 

/* first guess a format, then read `var' through `name',  macro `missing' 
 * is called if the entry is not found in the configuration file   */
#define CFG_GET_(var, name, tp, def, missing) {                       \
  char *fmt_ = CFG_TP2FMT(tp);                                        \
  if (fmt_[1] == '%') {  /* unable to determine format */             \
    fprintf(stderr, "cannot determine format for %s, "                \
        "file: %s, line: %d\n", #tp, __FILE__, __LINE__);             \
    goto ERR;                                                         \
  } else if (sizeof(var) != sizeof(tp)) {                             \
    fprintf(stderr, "var. %s does not seem to be a %s, "              \
        "file: %s, line: %d\n", #var, #tp, __FILE__, __LINE__);       \
    goto ERR;                                                         \
  } else if (0 != cfgget(cfg, &(var), name, fmt_)) {                  \
    missing(var, name, tp, def)                                       \
  }  }

/* variable to structure member */
#define CFG_STMBR(var)  (CFG_STPTR->var)
/* variable to key in configuration file */
#define CFG_STKEY(var)   CFG_STPFX #var

/* get structure member var, error if missing */
#define CFG_SMGETMUST(var, tp)                                        \
  CFG_GET_(CFG_STMBR(var), CFG_STKEY(var), tp, 0, CFG_MISERROR)

/* manually supply an other name */
#define CFG_SMGETMUST1(var, name, tp)                                 \
  CFG_GET_(CFG_STMBR(var), CFG_STKEY(name), tp, 0, CFG_MISERROR)

/* get structure member var, assume default if missing 
 * for string, def should be NULL or ssdup(...) */
#define CFG_SMGETDFLT(var, tp, def)                                   \
  CFG_STMBR(var) = def;                                               \
  CFG_GET_(CFG_STMBR(var), CFG_STKEY(var), tp, def, CFG_MISDFLT)

/* get structure member if cond is true, or def is set */
#define CFG_SMGETOPTN(var, tp, def, cond)                             \
  CFG_STMBR(var) = def;                                               \
  if (cond)                                                           \
    CFG_GET_(CFG_STMBR(var), CFG_STKEY(var), tp, def, CFG_MISDFLT)
  

/* demand structure member if cond is true, or def is set */
#define CFG_SMGETOPTM(var, tp, def, cond)                             \
  CFG_STMBR(var) = def;                                               \
  if (cond)                                                           \
    CFG_GET_(CFG_STMBR(var), CFG_STKEY(var), tp, def, CFG_MISERROR)

#define CFG_TESTLT(var, ref)                                          \
  if ((var) >= (ref)) {                                               \
    fprintf(stderr, "invalid value %s >= %s, "                        \
        "file: %s, line: %d\n", #var, #ref, __FILE__, __LINE__);      \
    goto ERR;                                                         \
  }

#define CFG_TESTGT(var, ref)                                          \
  if ((var) <= (ref)) {                                               \
    fprintf(stderr, "invalid value %s <= %s, "                        \
        "file: %s, line: %d\n", #var, #ref, __FILE__, __LINE__);      \
    goto ERR;                                                         \
  }

#define CFG_TESTIN(var, min, max)                                     \
  if ((var) < min || (var) > max) {                                   \
    fprintf(stderr, "invalid value %s outof (%s, %s), "               \
      "file: %s, line: %d\n", #var, #min, #max, __FILE__, __LINE__);  \
    goto ERR;                                                         \
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


