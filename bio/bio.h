#include "endn.h"
#ifndef BIO_H__
#define BIO_H__

#include <stdio.h>

/* 
 * Helper macros for reading binary files with endianness
 * support.  However, sizeof(int) must remain the same
 * between system and file.
 *
 * To use these macros in a function:
 * 1. define the following variables in your function
 *   FILE *fp;
 *   int endn, err;
 *   (no need for to define `endn' or `err' in writing a file)
 *
 * 2. define a label ERR for error exit
 *
 * 3. in reading a file, use BIO_INITENDIAN to determine 
 *    the correct endianness
 * */

#ifndef BIO_ENDNDEF
#define BIO_ENDNDEF 1  /* big endian */
#endif

/* string for printing file name and line number */
#define BIO_FLFMT_ "file: %s, line: %d"

/* check type */
#define BIO_CHECKTP_(x, tp)                                           \
  if (sizeof(x) != sizeof(tp)) {                                      \
    fprintf(stderr, "%s is not %s\n", #x, #tp);                       \
    goto ERR;                                                         \
  }

/* initialize file endian state to variable 'endn'
 * endn = 1: a conversion is needed from file's endianess to system's
 * endn = 0: otherwise
 * read an int variable x, 
 * determine endian by comparing the value of x with ref
 * quit if neither endians makes x == ref */
#define BIO_INIT_ENDIAN(x, ref) {                                     \
  BIO_CHECKTP_(x, int)                                                \
  if ((endn = endn_rmatchi(&(x), ref, fp)) < 0) {                     \
    fprintf(stderr, "%s 0x%X cannot match %s 0x%X\n",                 \
      #x, (unsigned) x, #ref, (unsigned) ref);                        \
    goto ERR;                                                         \
  } }

/* read an array of size n, set err, fix endian */
#define BIO_RATOM_(arr, n)                                            \
  if ((n) > 0 && endn_fread(arr, sizeof(*(arr)), n, fp, endn) != n) { \
    fprintf(stderr, "error while reading %s, size %u, "               \
        BIO_FLFMT_ "\n", #arr, (unsigned) n, __FILE__, __LINE__);     \
    err = 1;                                                          \
  } else { err = 0; }

/* read an array, set error */
#define BIO_RNA_(arr, n, tp) BIO_CHECKTP_(*(arr), tp) BIO_RATOM_(arr, n)
/* read a single variable x of type tp, set err if error occurs */
#define BIO_R1A_(x, tp) BIO_RNA_(&(x), 1, tp)

/* goto ERR if error occurs during reading */
#define BIO_RNB_(arr, n, tp) { BIO_RNA_(arr, n, tp); if (err) goto ERR; }
#define BIO_R1B_(x, tp) { BIO_R1A_(x, tp); if (err) goto ERR; }

/* most common: int and double cases */
#define BIO_RI_ERR(x)     BIO_R1A_(x, int)
#define BIO_RI(x)         BIO_R1B_(x, int)
#define BIO_RIARR(x, n)   BIO_RNB_(x, n, int)

#define BIO_RD_ERR(x)     BIO_R1A_(x, double)
#define BIO_RD(x)         BIO_R1B_(x, double)
#define BIO_RDARR(x, n)   BIO_RNB_(x, n, double)

/* match a temperory int x with the reference var */
#define BIO_MI(x, var)                                                \
  if ((x) != (var)) {                                                 \
    fprintf(stderr, "%s mismatch, expect: %d, read: %d "              \
        BIO_FLFMT_ "\n", #var, (int) var, x, __FILE__, __LINE__);     \
    goto ERR; }

/* match a temperory double x with the reference var */
#define BIO_MD(x, var, eps)                                           \
  if (fabs((x) - (var)) > eps) {                                      \
    fprintf(stderr, "%s mismatch, expect: %g, read: %g "              \
        BIO_FLFMT_ "\n", #var, var, x, __FILE__, __LINE__);           \
    goto ERR; }

/* read an int to x, match it with xref */
#define BIO_RMI(x, xref)       BIO_RI(x); BIO_MI(x, xref)
/* read a double to x, match it with xref */
#define BIO_RMD(x, xref, eps)  BIO_RD(x); BIO_MD(x, xref, eps)

/* write an array of size n with endian being BIO_ENDNDEF
 * we do not set err, directly goto ERR */
#define BIO_WATOM_(arr, n)                                            \
  if ((n) > 0 &&                                                      \
      endn_fwrite(arr, sizeof(*(arr)), n, fp, BIO_ENDNDEF) != n) {    \
    fprintf(stderr, "error while reading %s, size %u, "               \
        BIO_FLFMT_ "\n", #arr, (unsigned) n, __FILE__, __LINE__);     \
    goto ERR;                                                         \
  }

/* write an array, go to ERR if error occurs */
#define BIO_WNB_(arr, n, tp) BIO_CHECKTP_(*(arr), tp) BIO_WATOM_(arr, n)
/* write a single variable, go to ERR if error occurs */
#define BIO_W1B_(x, tp) BIO_WNB_(&(x), 1, tp)

#define BIO_WI(x)           BIO_W1B_(x, int)
#define BIO_WIARR(x, n)     BIO_WNB_(x, n, int)
#define BIO_WD(x)           BIO_W1B_(x, double)
#define BIO_WDARR(x, n)     BIO_WNB_(x, n, double)

#endif

