/*
  commonly-used routines
  Copyright (c) 2006-2013 Cheng Zhang

  This material is provided "as is", with absolutely no warranty expressed
  or implied. Any use is at your own risk.
  Permission to use or copy this software for any purpose is hereby granted
  without fee. Permission to modify the code and to distribute modified
  code is also granted without any restrictions.

  Usage:

  1.  It is designed for quick programming.
      For simple use, include this file and all functions will be available.
      But there might be many compiler warnings for unused functions.

  2.  You can include this file multiple times in a single file.

  3.  Function are static by default. To export functions,
      e.g., to make it easier to debug, or to avoid warnings of unused functions,
      define ZCOM_XFUNCS before the first inclusion.

  4.  To hand-pick specific set of modules, e.g.,
        #define ZCOM_PICK
        #define ZCOM_RNG
        #define ZCOM_ARGOPT
      before including this file, so other modules are skipped.

  5.  If the compiler supports keywords inline and restrict, write
        #define INLINE inline
        #define RESTRICT restrict
      before including this file. Otherwise the two keywords are guessed
      according to the compiler.

  6.  Define HAVEVAM if the compiler supports variable-argument macros.

  7.  The def module defines `real' as a double, to override it, write
        typedef float real;
        #define HAVEREAL 1
      before including this file (or equivalently define HAVE_REAL)
*/

/* ZCOM_PICK or ZCOM_NONE is used include only subset of modules
 * 1. to reduce the number of warnings for unused functions
 * 2. to reduce the compiling time
 * 3. to avoid potential name conflicts
 * By default, ZCOM_PICK is undefined, so everything is used. */
#ifdef ZCOM_NONE  /* equivalent to ZCOM_PICK */
#define ZCOM_PICK
#endif

#ifndef ZCOM_PICK
#endif

/* build dependencies */

/* manage storage class: static is the safer choice
   to avoid naming conflict.  Example:
   both m.c and n.c include this file,
   m.c --> m.o, n.c --> n.o, m.o+n.o --> a.out
   static is the only way to avoid naming conflict in this case.

   In case that this file is included multiple times,
   ZCOM_XFUNCS should be defined before the first inclusion,
   otherwise it won't be effective in deciding storage class. */
#ifndef STRCLS
  #ifndef ZCOM_XFUNCS
    #define STRCLS static
  #else
    #define STRCLS
  #endif
#endif

/* inline keyword */
#ifndef INLINE
  #if defined(__GNUC__) || defined(__xlC__)
    #define INLINE STRCLS __inline__
  #elif defined(__INTEL_COMPILER)
    #define INLINE STRCLS __inline
  #elif defined(_MSC_VER) || defined(__BORLANDC__)
    #define INLINE __inline STRCLS
  #elif defined(__STDC_VERSION__) && (STDC_VERSION__ >= 199901L)
    #define INLINE STRCLS inline
  #else
    #define INLINE STRCLS
  #endif
#endif

/* restrict keyword */
#ifndef RESTRICT
  #if (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__xlC__) \
      || (defined(_MSC_VER) && _MSC_VER >= 1400))
    #define RESTRICT __restrict
  #elif defined(__STDC_VERSION__) && (STDC_VERSION__ >= 199901L)
    #define RESTRICT restrict
  #else
    #define RESTRICT
  #endif
#endif

/* macros with variable-length arguments */
#ifndef HAVEVAM
  #if (  (defined(__GNUC__) && (__GNUC__ >= 3))   \
      || (defined(__xlC__)  && (__xlC__ >= 0x0700)) \
      || (defined(_MSC_VER) && (_MSC_VER >= 1400)) )
    #define HAVEVAM 1
  #endif
#endif

#ifdef __INTEL_COMPILER
  #pragma warning push
  #pragma warning(disable:981) /* unspecified order of operands */
  #pragma warning(disable:177) /* unreferenced function */
  #pragma warning(disable:161) /* unrecognized #pragma, for omp */
#elif defined(__GNUC__) && (__GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
  /* Note: intel compiler also defines __GNUC__ */
  #if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
    /* diagonistic push and pop are added in GCC 4.6 */
    #pragma GCC diagnostic push
  #endif
  #pragma GCC diagnostic ignored "-Wunknown-pragmas"
  #pragma GCC diagnostic ignored "-Wvariadic-macros"
#endif

#ifdef _MSC_VER
  #pragma warning(push)
  #pragma warning(disable:4127) /* conditional expression constant */
  #pragma warning(disable:4505) /* unreferenced function */
  #pragma warning(disable:4514) /* unreferenced inline */
  #pragma warning(disable:4710) /* not inlined */
  /* suppress CRT safety warnings:
   * no need for "_s" versions of CRT functions
   * The macros below will not work if standard C library files
   * like <stdio.h> and <time.h> are included before zcom.h */
  #ifndef _CRT_SECURE_NO_DEPRECATE
  #define _CRT_SECURE_NO_DEPRECATE
  #endif
  #ifndef _CRT_SECURE_NO_WARNINGS
  #define _CRT_SECURE_NO_WARNINGS
  #endif
  #include <stdio.h>
#endif

#ifdef __BORLANDC__
  #pragma warning(push)
  #pragma warn -8027 /* not expanded inlined */
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* In addition to ZCOM_ABC, we have to define another macro ZCOM_ABC__
 * to allow multiple inclusions. */


/* finishing up */
#ifdef __INTEL_COMPILER
  #pragma warning pop
#elif defined(__GNUC__) && ( __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6))
  /* diagonistic push and pop are added in GCC 4.6 */
  #pragma GCC diagnostic pop
#endif

#ifdef _MSC_VER
  #pragma warning(pop)
#endif
