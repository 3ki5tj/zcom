/*
  A set of common (but not so optimized) routines.

  Copyright (c) 2006-2010 Cheng Zhang

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Usage:

  1.  This file is mainly intended to be included directly,
      not to be compiled independently.

  2.  It is designed to facilitate simple programming.
      For simple use, just include it with no extra pre-definition.
      All functions will be available.
      But there might be compiler warnings for unused functions.

  3.  You should be able to include this file multiple times in a single file
      without a problem (otherwise a bug!).

  4.  Every function is static by default. If you want to export functions,
      e.g., to make it easier to debug, or to avoid warning of unused functions,
      define ZCOM_XFUNCS before the first inclusion.

  5.  The file is a combo of separate modules, to include just one module,
      define ZCOM_PICK and ZCOM_xxx before inclusion, e.g.
        #define ZCOM_PICK
        #define ZCOM_RNG
        #include "name_of_this_file"

      One can probably save some compiling time and avoid some warnings
      in this way. It also makes more sense.

  6.  Define ZCHAVEVAM if the compiler supports variable-argument macros.

  7.  A few type-names defined depending the modules used,
      e.g. cfgdata_t if ZCOM_CFG is defined.

  -------------------------------------------------------------------------------

  8.  By very careful about the RV3/RV2 module, because it will typedef `real',
      and by default to double.
      If you want to change it define ZCREAL as float.
      If you already defined real somewhere else (as float or double),
      you must
        #define ZCHAVEREAL
      before including this file.
      Further, to be compatible with GROMACS, you can define HAVE_REAL to achive
      the same effect.
      The same applies to the two-dimensional version RV2

  9.  When using RV2/RV3 module, you can further use
        #define ZCOM_RVDIM_BINDING 2
      before inclusion, so the rv_xxxx is mapped to rv2_xxxx.
      The same applies to the 3d version.

      However, these mapping
      + are not setup by default.
      + only applies to routines common to 2d and 3d.
*/


/* A few macros to allow user to selectively include functions.
   It helps
   1. to reduce the # of warnings for unused functions
   2. to accelerate the compiling
   3. avoid multiple inclusions
   By default (ZCOM_PICK not defined), we use everything old
 */
#ifdef ZCOM_NONE  // equivalent to ZCOM_PICK
#define ZCOM_PICK
#endif

#ifndef ZCOM_PICK
  #ifndef ZCOM_DIE
  #define ZCOM_DIE
  #endif
  #ifndef ZCOM_SS
  #define ZCOM_SS
  #endif
  #ifndef ZCOM_RNG
  #define ZCOM_RNG
  #endif
  #ifndef ZCOM_CFG
  #define ZCOM_CFG
  #endif
  #ifndef ZCOM_TRACE
  #define ZCOM_TRACE
  #endif
  #ifndef ZCOM_DIST
  #define ZCOM_DIST
  #endif
  #ifndef ZCOM_LOG
  #define ZCOM_LOG
  #endif
  #ifndef ZCOM_STR
  #define ZCOM_STR
  #endif
  #ifndef ZCOM_LU
  #define ZCOM_LU
  #endif
  #ifndef ZCOM_ZT
  #define ZCOM_ZT
  #endif
  #ifndef ZCOM_RV3
  #define ZCOM_RV3
  #endif
  #ifndef ZCOM_RV2
  #define ZCOM_RV2
  #endif
  #ifndef ZCOM_DIHCALC
  #define ZCOM_DIHCALC
  #endif
  #ifndef ZCOM_ENDN
  #define ZCOM_ENDN
  #endif
  #ifndef ZCOM_XM
  #define ZCOM_XM
  #endif
#endif

/* build dependencies */
#if (defined(ZCOM_RNG)  || defined(ZCOM_TRACE) || defined(ZCOM_CFG) || \
     defined(ZCOM_HIST) || defined(ZCOM_LOG)   || defined(ZCOM_ZT) )
  #define ZCOM_SS  /* needs file name support */
#endif

#ifdef ZCOM_XM
  #define ZCOM_SS
  #define ZCOM_CFG
  #define ZCOM_DIE
#endif

#ifdef ZCOM_CFG
  #define ZCOM_DIE
  #define ZCOM_SS
#endif

#ifdef ZCOM_SS
  #define ZCOM_DIE
#endif

#ifdef ZCOM_DIHCALC
  #define ZCOM_RV3
#endif

#ifdef ZCOM_BIO
  #define ZCOM_ENDN
#endif

/* manage storage class: static is still the safer choice
   to avoid naming conclict.  Example:
   both m.c and n.c include this file,
   m.c --> m.o, n.c --> n.o, m.o+n.o --> a.out
   static is the only way to avoid naming conflict in this case.

   In case that this file is included multiple times,
   ZCOM_XFUNCS should be defined before the first inclusion,
   otherwise it won't be effective in deciding storage class.
 */
#ifndef ZCOM_XFUNCS
  #ifndef ZCSTRCLS
    #define ZCSTRCLS static
  #endif
  #ifndef ZCINLINE
    #define ZCINLINE static __inline
  #endif
#else
  #ifndef ZCSTRCLS
    #define ZCSTRCLS
  #endif
  #ifndef ZCINLINE
    /* should be `inline' according to C99,
       but gcc, vc++, intel support `__inline' by default */
    #define ZCINLINE __inline
  #endif
#endif

/* try to use restrict if possible */
#ifndef ZCRESTRICT
  #if (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__xlC__))
    /* should be `restrict' according to C99,
       but compilers, e.g., gcc and icc, support `__restrict' by default */
    #define ZCRESTRICT __restrict
  #else
    #define ZCRESTRICT
  #endif
#endif

/* newer compilers should support macros with variable-length arguments */
#ifndef ZCHAVEVAM
  #if ( (defined(__GNUC__)&&(__GNUC__>=3)) || (defined(__xlC__)&&(__xlC__>=0x0700)) || (defined(_MSC_VER)&&(_MSC_VER>=1400)) )
    #define ZCHAVEVAM 1
  #endif
#endif

/* In addition to ZCOM_ABC, we have to define another macro ZCOM_ABC__
 * in order to avoid multiple inclusion a single ZCOM_ABC__ won't do,
 * because different module-set may be selected */
#ifdef  ZCOM_DIE
#ifndef ZCOM_DIE__
#define ZCOM_DIE__

#endif /* ZCOM_DIE__ */
#endif /* ZCOM_DIE */


#ifdef  ZCOM_SS
#ifndef ZCOM_SS__
#define ZCOM_SS__

#endif /* ZCOM_SS__ */
#endif /* ZCOM_SS */


#ifdef  ZCOM_RNG
#ifndef ZCOM_RNG__
#define ZCOM_RNG__

#endif /* ZCOM_RNG__ */
#endif /* ZCOM_RNG */

#ifdef  ZCOM_CFG
#ifndef ZCOM_CFG__
#define ZCOM_CFG__

#endif /* ZCOM_CFG__ */
#endif /* ZCOM_CFG */

#ifdef  ZCOM_TRACE
#ifndef ZCOM_TRACE__
#define ZCOM_TRACE__
#endif /* ZCOM_TRACE__ */
#endif /* ZCOM_TRACE */


#ifdef  ZCOM_DIST
#ifndef ZCOM_DIST__
#define ZCOM_DIST__

#endif /* ZCOM_DIST__ */
#endif /* ZCOM_DIST */


#ifdef  ZCOM_LOG
#ifndef ZCOM_LOG__
#define ZCOM_LOG__

#endif /* ZCOM_LOG__ */
#endif /* ZCOM_LOG */


#ifdef  ZCOM_LU
#ifndef ZCOM_LU__
#define ZCOM_LU__

#endif /* ZCOM_LU__ */
#endif /* ZCOM_LU */

#ifdef  ZCOM_STR
#ifndef ZCOM_STR__
#define ZCOM_STR__

#include <ctype.h>
#include <stdio.h>

#define ZCOM_STRCMP  0x0001
#define ZCOM_STRCPY  0x0002
#define ZCOM_STRCAT  0x0004

#define ZCOM_DOCASE  0x0100
#define ZCOM_UPPER   0x0200

/* copy the string and convert it to upper/lower case */
#define strconv2upper(s,t,size_s) zcom_strconv(s,t,size_s,ZCOM_STRCPY|ZCOM_DOCASE|ZCOM_UPPER)
#define strconv2lower(s,t,size_s) zcom_strconv(s,t,size_s,ZCOM_STRCPY|ZCOM_DOCASE)
#define strcpy_safe(s,t,size_s)   zcom_strconv(s,t,size_s,ZCOM_STRCPY)
/* concatenate strings, the last parameter is the buffer size of s,
 * unlike strncat(), in which it's the number of characters from *t* to be copied.  */
#define strcat_safe(s,t,size_s)   zcom_strconv(s,t,size_s,ZCOM_STRCAT)
/* compare strings, ignore case differences, stricmp */
#define strcmp_igncase(s,t)       zcom_strconv((char *)s,t,0,ZCOM_STRCMP|ZCOM_DOCASE)

/* A combination of string comparison and conversion.
 * It can compare strings with or without cases,
 * and it can convert strings to different cases.
 *
 * To balance the string comparison and copying,
 * we make the return value int, instead of char *.
 *
 * In case of string-copying, it returns
 * 0 for success,
 * 1 for lack of space, or NULL pointers.
 * Unlike the case of strncpy(), s is always null-terminated.
 * */
int zcom_strconv(char *s, const char *t, size_t size_s, unsigned flags)
{
  size_t i, j;
  int cs, ct, docase, doupper;

  docase = (flags&ZCOM_DOCASE); /* do case conversion */
  doupper = (flags&ZCOM_UPPER);

  if (flags&ZCOM_STRCMP) { /* comparison, size_s is ignored */
    if (s == NULL || t == NULL) return 0;
    for (i = 0; ; i++) {
      cs = s[i];
      ct = t[i];
      if (docase) {
        cs = tolower((unsigned char)cs);
        ct = tolower((unsigned char)ct);
      }
      if (cs == 0 || ct == 0 || cs != ct) break;
    }
    return cs-ct;
  } else if (flags & (ZCOM_STRCPY|ZCOM_STRCAT)) { /* copying and */
    if (size_s == 0 || s == NULL || t == NULL) return 1;
    /* t[size_s-1] should be '\0' */
    i = 0;
    if (flags&ZCOM_STRCAT) while(s[i]) i++;
    for (j = 0; i < size_s-1; i++, j++) {
      ct = t[j];
      if (docase && (ct != 0)) {
        ct = (unsigned char)(char)ct;
        if (doupper) {
          cs = toupper(ct);
        } else {
          cs = tolower(ct);
        }
      } else {
        cs = ct;
      }
      s[i] = (char)cs;
      if (ct == 0) break;
    }
    if (i == size_s-1) s[i] = '\0';
    return (t[j] != '\0');
  } else { /* unknown flag */
    fprintf(stderr, "zcom_strconv: invalid flags=%#x.\n", flags);
    return 0;
  }
}

#endif /* ZCOM_STR__ */
#endif /* ZCOM_STR */


#ifdef  ZCOM_ZT
#ifndef ZCOM_ZT__
#define ZCOM_ZT__

#endif /* ZCOM_ZT__ */
#endif /* ZCOM_ZT */

#ifdef ZCOM_RV3
#ifndef ZCOM_RV3__
#define ZCOM_RV3__

#ifdef HAVE_REAL
  #ifndef ZCHAVEREAL
  #define ZCHAVEREAL HAVE_REAL
  #endif
#endif

/* Usage: if you already typedef'ed real, then define ZCHAVEREAL
   before including this file to avoid a conflict.
   Otherwise define ZCREAL as float or double */
#ifndef ZCHAVEREAL
  #define ZCHAVEREAL 1
  #ifndef ZCREAL
    #define ZCREAL double
  #endif
  typedef ZCREAL real;
#endif

#include <math.h>

/* due to that pointer may overlap with each other,
 * be careful when using the const modifier */

ZCINLINE void rv3_zero(real *x) { x[0] = 0.0f; x[1] = 0.0f; x[2] = 0.0f; }
ZCINLINE void rv3_copy(real *x, const real *src) { x[0] = src[0]; x[1] = src[1]; x[2] = src[2]; }

ZCINLINE real rv3_sqr (const real *x) { return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; }
ZCINLINE real rv3_norm(const real *x) { return (real)sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

ZCINLINE real *rv3_normalize(real *x)
{
  real r = rv3_norm(x);
  if (r > 0.0) {
    r = 1.0f/r;
    x[0] *= r;
    x[1] *= r;
    x[2] *= r;
  }
  return x;
}

/* if x == y, try to use sqr */
ZCINLINE real rv3_dot(const real *x, const real *y)
{
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

ZCINLINE real *rv3_cross(real *ZCRESTRICT z, const real *x, const real *y)
{
  z[0] = x[1]*y[2]-x[2]*y[1];
  z[1] = x[2]*y[0]-x[0]*y[2];
  z[2] = x[0]*y[1]-x[1]*y[0];
  return z;
}

ZCINLINE real *rv3_neg(real *x)
{
  x[0] -= x[0];
  x[1] -= x[1];
  x[2] -= x[2];
  return x;
}

ZCINLINE real *rv3_neg2(real *nx, const real *x)
{
  nx[0] = -x[0];
  nx[1] = -x[1];
  nx[2] = -x[2];
  return nx;
}

ZCINLINE real *rv3_inc(real *x, const real *dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  x[2] += dx[2];
  return x;
}
ZCINLINE real *rv3_dec(real *x, const real *dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  x[2] -= dx[2];
  return x;
}
ZCINLINE real *rv3_sinc(real *x, real *dx, real s)
{
  x[0] += s*dx[0];
  x[1] += s*dx[1];
  x[2] += s*dx[2];
  return x;
}
ZCINLINE real *rv3_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return x;
}

/* if y == x, just use smul */
ZCINLINE real *rv3_smul2(real * ZCRESTRICT y, const real *x, real s)
{
  y[0] = x[0]*s;
  y[1] = x[1]*s;
  y[2] = x[2]*s;
  return y;
}

/* for in-place difference use rv3_dec */
ZCINLINE real *rv3_diff(real * ZCRESTRICT diff, const real *a, const real *b)
{
  diff[0] = a[0]-b[0];
  diff[1] = a[1]-b[1];
  diff[2] = a[2]-b[2];
  return diff;
}

/* sum = a+b, for in-place addition use rv3_inc */
ZCINLINE real *rv3_sum2(real * ZCRESTRICT sum, const real *a, const real *b)
{
  sum[0] = a[0]+b[0];
  sum[1] = a[1]+b[1];
  sum[2] = a[2]+b[2];
  return sum;
}

/* sum = -a-b */
ZCINLINE real *rv3_nsum2(real *sum, const real *a, const real *b)
{
  sum[0] = -a[0]-b[0];
  sum[1] = -a[1]-b[1];
  sum[2] = -a[2]-b[2];
  return sum;
}

ZCINLINE real *rv3_lincomb2(real *sum, const real *a, const real *b, real s1, real s2)
{
  sum[0] = a[0]*s1+b[0]*s2;
  sum[1] = a[1]*s1+b[1]*s2;
  sum[2] = a[2]*s1+b[2]*s2;
  return sum;
}
#endif /* ZCOM_RV3__ */
#endif /* ZCOM_RV3 */

#ifdef  ZCOM_DIHCALC
#ifndef ZCOM_DIHCALC__
#define ZCOM_DIHCALC__

#endif /* ZCOM_CALCDIH__ */
#endif /* ZCOM_CALCDIH */

#ifdef ZCOM_RV2
#ifndef ZCOM_RV2__
#define ZCOM_RV2__

#ifdef HAVE_REAL
  #ifndef ZCHAVEREAL
  #define ZCHAVEREAL HAVE_REAL
  #endif
#endif

/* Usage: if you already typedef'ed real, then define ZCHAVEREAL
   before including this file to avoid a conflict.
   Otherwise define ZCREAL as float or double */
#ifndef ZCHAVEREAL
  #define ZCHAVEREAL 1
  #ifndef ZCREAL
    #define ZCREAL double
  #endif
  typedef ZCREAL real;
#endif

#include <math.h>

/* due to that pointer may overlap with each other,
 * be careful when using the const modifier */

ZCINLINE void rv2_zero(real *x) { x[0] = 0; x[1] = 0; }
ZCINLINE void rv2_copy(real *x, const real *src) { x[0] = src[0]; x[1] = src[1]; }

ZCINLINE real rv2_sqr(const real *x) { return x[0]*x[0]+x[1]*x[1]; }
ZCINLINE real rv2_norm(const real *x) { return (real)sqrt(x[0]*x[0]+x[1]*x[1]); }

ZCINLINE real *rv2_normalize(real *x)
{
  real r = rv2_norm(x);
  if (r > 0.0) {
    r = 1.0/r;
    x[0] *= r;
    x[1] *= r;
  }
  return x;
}

ZCINLINE real rv2_dot(const real *x, const real *y) { return x[0]*y[0]+x[1]*y[1]; }

ZCINLINE real rv2_cross(const real *x, const real *y)
{
  return x[0]*y[1]-x[1]*y[0];
}

ZCINLINE real *rv2_neg(real *x)
{
  x[0] -= x[0];
  x[1] -= x[1];
  return x;
}

ZCINLINE real *rv2_neg2(real *nx, const real *x)
{
  nx[0] = -x[0];
  nx[1] = -x[1];
  return nx;
}

ZCINLINE real *rv2_inc(real *x, const real *dx)
{
  x[0] += dx[0];
  x[1] += dx[1];
  return x;
}
ZCINLINE real *rv2_dec(real *x, const real *dx)
{
  x[0] -= dx[0];
  x[1] -= dx[1];
  return x;
}
ZCINLINE real *rv2_sinc(real *x, real *dx, real s)
{
  x[0] += s*dx[0];
  x[1] += s*dx[1];
  return x;
}
ZCINLINE real *rv2_smul(real *x, real s)
{
  x[0] *= s;
  x[1] *= s;
  return x;
}
ZCINLINE real *rv2_smul2(real *y, const real *x, real s)
{
  y[0] = x[0]*s;
  y[1] = x[1]*s;
  return y;
}

/* for in-place difference use rv3_dec */
ZCINLINE real *rv2_diff(real *diff, const real *a, const real *b)
{
  diff[0] = a[0]-b[0];
  diff[1] = a[1]-b[1];
  return diff;
}

/* sum = a+b, for in-place addition use rv3_inc */
ZCINLINE real *rv2_sum2(real *sum, const real *a, const real *b)
{
  sum[0] = a[0]+b[0];
  sum[1] = a[1]+b[1];
  return sum;
}

/* sum = -a-b */
ZCINLINE real *rv2_nsum2(real *sum, const real *a, const real *b)
{
  sum[0] = -a[0]-b[0];
  sum[1] = -a[1]-b[1];
  return sum;
}

ZCINLINE real *rv2_lincomb2(real *sum, const real *a, const real *b, real s1, real s2)
{
  sum[0] = a[0]*s1+b[0]*s2;
  sum[1] = a[1]*s1+b[1]*s2;
  return sum;
}

#endif /* ZCOM_RV2__ */
#endif /* ZCOM_RV2 */


#if (defined(ZCOM_RVDIM_BINDING) && !defined(ZCOM_RVDIM_BINDING_DEFINED))
/* these binding macros reside outside both RV2 and RV3
 * to avoid undesirable shielding during multiple inclusion, e.g.,
 * in the first inclusion ZCOM_RVDIM_BINDING is not defined,
 * but it is in the second one.
 */
  #define ZCOM_RVDIM_BINDING_DEFINED 1
  #if (ZCOM_RVDIM_BINDING==2)
    #define rv_zero       rv2_zero
    #define rv_copy       rv2_copy
    #define rv_sqr        rv2_sqr
    #define rv_norm       rv2_norm
    #define rv_normalize  rv2_normalize
    #define rv_dot        rv2_dot
    #define rv_cross      rv2_cross
    #define rv_neg        rv2_neg
    #define rv_neg2       rv2_neg2
    #define rv_inc        rv2_inc
    #define rv_dec        rv2_dec
    #define rv_sinc       rv2_sinc
    #define rv_smul       rv2_smul
    #define rv_smul2      rv2_smul2
    #define rv_diff       rv2_diff
    #define rv_sum2       rv2_sum2
    #define rv_nsum2      rv2_nsum2
    #define rv_lincomb2   rv2_lincomb2
  #else
    #define rv_zero       rv3_zero
    #define rv_copy       rv3_copy
    #define rv_sqr        rv3_sqr
    #define rv_norm       rv3_norm
    #define rv_normalize  rv3_normalize
    #define rv_dot        rv3_dot
    #define rv_cross      rv3_cross
    #define rv_neg        rv3_neg
    #define rv_neg2       rv3_neg2
    #define rv_inc        rv3_inc
    #define rv_dec        rv3_dec
    #define rv_sinc       rv3_sinc
    #define rv_smul       rv3_smul
    #define rv_smul2      rv3_smul2
    #define rv_diff       rv3_diff
    #define rv_sum2       rv3_sum2
    #define rv_nsum2      rv3_nsum2
    #define rv_lincomb2   rv3_lincomb2
  #endif
#endif

#ifdef  ZCOM_ENDN
#ifndef ZCOM_ENDN__
#define ZCOM_ENDN__

#endif /* ZCOM_ENDN__ */
#endif /* ZCOM_ENDN */

#ifdef  ZCOM_BIO
#ifndef ZCOM_BIO__
#define ZCOM_BIO__

#endif /* ZCOM_BIO__ */
#endif /* ZCOM_BIO */

#ifdef  ZCOM_XM
#ifndef ZCOM_XM__
#define ZCOM_XM__

#endif /* ZCOM_XM__ */
#endif /* ZCOM_XM */

