/*
  common routines

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
        #include "this_file"

  6.  Define ZCHAVEVAM if the compiler supports variable-argument macros.

  7.  A few type-names defined depending the modules used,
      e.g. cfgdata_t if ZCOM_CFG is defined.

  -------------------------------------------------------------------------------

  8.  By careful about the RV3/RV2 module, because it will typedef `real'
      to double, to override: 
        typedef float/double real;
        #define ZCHAVEREAL 1
      before including this file (or equivalently HAVE_REAL)

  9.  When using RV2/RV3 module, you can further use
        #define ZCOM_RVDIM_BINDING 2
      before inclusion, so the rv_xxxx is mapped to rv2_xxxx.
      The same applies to the 3d version.

      However, these mapping
      + are not setup by default.
      + only applies to routines common to 2d and 3d.
*/

/* macros for selectively including functions, advantages: 
 * 1. reduces the # of warnings for unused functions
 * 2. accelerates the compiling
 * 3. avoids multiple inclusions
 * By default (ZCOM_PICK not defined), everything is used */
#ifdef ZCOM_NONE  /* equivalent to ZCOM_PICK */
#define ZCOM_PICK
#endif

#ifndef ZCOM_PICK
#endif

/* build dependencies */

/* manage storage class: static is the safer choice
   to avoid naming conclict.  Example:
   both m.c and n.c include this file,
   m.c --> m.o, n.c --> n.o, m.o+n.o --> a.out
   static is the only way to avoid naming conflict in this case.

   In case that this file is included multiple times,
   ZCOM_XFUNCS should be defined before the first inclusion,
   otherwise it won't be effective in deciding storage class. */
#ifndef ZCSTRCLS
  #ifndef ZCOM_XFUNCS
    #define ZCSTRCLS static
  #else
    #define ZCSTRCLS
  #endif
#endif

/* inline keyword */
#ifndef ZCINLINE
  #if defined(__GNUC__) || defined(__xlC__)
    #define ZCINLINE ZCSTRCLS __inline__
  #elif defined(_MSC_VER) || defined(__BORLANDC__)
    #define ZCINLINE __inline ZCSTRCLS
  #elif defined(__STDC_VERSION__) && (STDC_VERSION__ >= 199901L)
    #define ZCINLINE ZCSTRCLS inline
  #else
    #define ZCINLINE ZCSTRCLS
  #endif
#endif

/* restrict keyword */
#ifndef ZCRESTRICT
  #if (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__xlC__))
    #define ZCRESTRICT __restrict
  #elif defined(__STDC_VERSION__) && (STDC_VERSION__ >= 199901L)
    #define ZCRESTRICT restrict
  #else
    #define ZCRESTRICT
  #endif
#endif

/* macros with variable-length arguments */
#ifndef ZCHAVEVAM
  #if (  (defined(__GNUC__) && (__GNUC__ >= 3))   \
      || (defined(__xlC__)  && (__xlC__ >= 0x0700)) \
      || (defined(_MSC_VER) && (_MSC_VER >= 1400)) ) 
    #define ZCHAVEVAM 1
  #endif
#endif

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
    #define rv_ncopy      rv2_ncopy
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
    #define rv_dist       rv2_dist
    #define rv_dist2      rv2_dist2
    #define rv_add        rv2_add
    #define rv_nadd       rv2_nadd
    #define rv_lincomb2   rv2_lincomb2
  #else
    #define rv_zero       rv3_zero
    #define rv_copy       rv3_copy
    #define rv_ncopy      rv3_ncopy
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
    #define rv_dist       rv3_dist
    #define rv_dist2      rv3_dist2
    #define rv_add        rv3_add
    #define rv_nadd       rv3_nadd
    #define rv_lincomb2   rv3_lincomb2
  #endif
#endif

/* In addition to ZCOM_ABC, we have to define another macro ZCOM_ABC__
 * in order to avoid multiple inclusion a single ZCOM_ABC__ won't do,
 * because different module-set may be selected */

