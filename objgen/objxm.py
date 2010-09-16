# content of helper macros
helper_macros = r'''
#ifndef X_MACROS__
#define X_MACROS__

#include <stdio.h>

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* return a format string from a type
 * undetermined format goto "%%"  */
#define XM_TP2FMT_(tp) \
    (strcmp(#tp, "int")       == 0 ? "%d"  : \
     strcmp(#tp, "unsigned")  == 0 ? "%u"  : \
     strcmp(#tp, "float")     == 0 ? "%f"  : \
     strcmp(#tp, "double")    == 0 ? "%lf" : \
     strcmp(#tp, "char *")    == 0 ? "%s"  : \
     strcmp(#tp, "null")      == 0 ? ""    : "%%")

#define XM_PRINT_FILE_LINE_() \
    fprintf(stderr, "file: %s, line: %d\n", __FILE__, __LINE__)

#ifndef XM_FATAL_ACTION_  /* for programmer's mistake */
#define XM_FATAL_ACTION_  { XM_PRINT_FILE_LINE_(); exit(1); }
#endif

/* print a message leading by `msg', and execute `action' 
 * low level error handling */
#define XM_ERRMSG_(var, key, tp, def, desc, msg, action) {            \
    fprintf(stderr, "%s: var: %s, key: %s, type: %s, def: %s\n",      \
        msg, #var, key, #tp, #def);                                   \
    /* if (desc[0] != '\0') fprintf(stderr, "desc: %s, ", desc); */   \
    /* fprintf(stderr, "def: %s\n\n", #def); */                       \
    /* XM_PRINT_FILE_LINE_() */                                       \
    action; }

/* first guess a format, then read `var' through `key' (char *),
 * if the variable is missing, mismsg is printed, with misact
 * is called if the entry is not found in the configuration file   
 * the whole process is skipped if null is passed to tp 
 * */
#define XM_CFGGET_(cfg, var, key, tp, def, desc, mismsg, misact) {    \
  char *fmt_ = XM_TP2FMT_(tp);                                        \
  if (fmt_[0] == '\0') {                                              \
    /* do nothing */ ;                                                \
  } else if (fmt_[1] == '%') {  /* unable to determine format */      \
    fprintf(stderr, "cannot determine format for %s\n", #tp);         \
    XM_FATAL_ACTION_; /* fatal: programmer's mistake */               \
  } else if (sizeof(var) != sizeof(tp)) {                             \
    fprintf(stderr, "var. %s is not of type %s\n", #var, #tp);        \
    XM_FATAL_ACTION_; /* fatal: programmer's mistake */               \
  } else if (0 != cfgget(cfg, &(var), key, fmt_)) {                   \
    XM_ERRMSG_(var, key, tp, def, desc, mismsg, misact);              \
  }  }

/* conditionally (t0) get var from configuration file, 
 * then check if condition t1 is met, if not, misact is executed
 * neither t0 nor t1 should be empty, but can be 1 or 0;
 * finally, expression `eval', which can be empty, is evaluated;
 * a default value 'def' is always assigned at the beginning for safety,
 * the design is based on the concern that if anything fails
 * an uninitialized variable can be dangerous
 * */
#define XM_CFGGETC_(cfg, var, key, tp, def, t0, t1, eval, desc, mismsg, misact) \
  var = def;                                                                    \
  if (XM_TEST_(t0)) {                                                           \
    XM_CFGGET_(cfg, var, key, tp, def, desc, mismsg, misact)                    \
    XM_TESTERR_(t1, misact)                                                     \
    eval;                                                                       \
  }

/* test if expression t is true
 * t should *not* be empty 
 * t can be 1; since we tested #t[0], compilers usually don't complain
 *   the condition is always true; a further safer choice is probably
 *   strcmp(#t, "1") == 0
 * */
#define XM_TEST_(t) ((#t[0] == '1' && #t[1] == '\0') ? 1 : (t))

/* if condition 'cond' is not true, execute `erract' */
#define XM_TESTERR_(t, erract)                                    \
  if ( !XM_TEST_(t) )  {                                          \
    fprintf(stderr, "failed cond: %s\n", #t);                     \
    XM_PRINT_FILE_LINE_();                                        \
    erract;                                                       \
  }

/* initialize a static array arr with def
 * `def' can be a number, or an expression using the index ia_  
 * `assign' can be empty for simple array initialization, but it can 
 * also be a function call/statement that initializes the whole array */
#define XM_INITARR_(arr, cnt, def, assign)               \
   if (#assign[0] == '\0') {                             \
    size_t ia_;                                          \
    for (ia_ = 0; ia_ < (cnt); ia_++) arr[ia_] = (def);  \
   } else { assign; }

/* allocate array, and initialize it with function def
 * tp should be a pointer type, like int *, or double * 
 * always allocate cnt+1 elements for safety, 
 * but only 0 to cnt-1 are initialized
 * */
#define XM_DARR_(arr, tp, def, cnt, assign, desc, errmsg, erract)           \
  if (sizeof(arr[0]) != sizeof(*(tp)arr)) {                                 \
    fprintf(stderr, "dynamic array %s is not of type %s\n", #arr, #tp);     \
    XM_ERRMSG_(arr, #arr, tp, def, desc, "type error", XM_FATAL_ACTION_);   \
  } else if ((arr = calloc( (cnt+1), sizeof(arr[0]) )) == NULL) {           \
    XM_ERRMSG_(arr, #arr, tp, def, desc, errmsg, erract);                   \
  } else {                                                                  \
    XM_INITARR_(arr, cnt, def, assign)                                      \
  }

/* initialize a static array */
#define XM_SARR_(arr, tp, def, cnt, assign, desc)                           \
  if (sizeof(arr[0]) != sizeof(tp)) {                                       \
    fprintf(stderr, "static array %s is not of type %s\n", #arr, #tp);      \
    XM_ERRMSG_(arr, #arr, tp, def, desc, "type error", XM_FATAL_ACTION_);   \
  } else {                                                                  \
    XM_INITARR_(arr, cnt, def, assign)                                      \
  }

#endif
'''


