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
  #ifndef ZCOM_ERROR
  #define ZCOM_ERROR
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
  #ifndef ZCOM_ENDIAN
  #define ZCOM_ENDIAN
  #endif
#endif

/* build dependencies */
#if (defined(ZCOM_RNG)  || defined(ZCOM_TRACE) || defined(ZCOM_CFG) || \
     defined(ZCOM_HIST) || defined(ZCOM_LOG)   || defined(ZCOM_ZT) )
  #define ZCOM_SS  /* needs file name support */
#endif

#ifdef ZCOM_SS
  #define ZCOM_ERROR
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

#ifdef  ZCOM_ERROR
#ifndef ZCOM_ERROR__
#define ZCOM_ERROR__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/* print an error message, then quit */
/* newer compilers should support macros with variable-length arguments */
#if ( (defined(__GNUC__)&&(__GNUC__>=3)) || (defined(__xlC__)&&(__xlC__>=0x0700)) || (defined(_MSC_VER)&&(_MSC_VER>=1400)) )
#define zcom_fatal(fmt, ...)  zcom_fatal_(__FILE__, __LINE__, fmt, ## __VA_ARGS__)
ZCSTRCLS void zcom_fatal_(const char *file, int lineno, const char *fmt, ...)
#else
ZCSTRCLS void zcom_fatal(const char *fmt, ...)
#endif
{
  va_list args;
  fprintf(stderr, "Fatal ");
  if(file != NULL) fprintf(stderr, "%s ", file);
  if(lineno > 0)   fprintf(stderr, "%d ", lineno);
  fprintf(stderr, "| ");
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  exit(1);
}

#endif /* ZCOM_ERROR__ */
#endif /* ZCOM_ERROR */


#ifdef  ZCOM_SS
#ifndef ZCOM_SS__
#define ZCOM_SS__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

enum { SSCAT=1, SSDELETE=2, SSSHRINK=3, SSSINGLE=0x1000 };

#define ssnew(n)       sscpycatx(NULL, NULL, (n),    0)
#define ssdup(t)       sscpycatx(NULL, (t),   0,     0)
#define sscpy(s, t)    sscpycatx(&(s), (t),   0,     0)
#define sscat(s, t)    sscpycatx(&(s), (t),   0, SSCAT)
#define ssdel(s)       ssmanage((s), SSDELETE|SSSINGLE)
#define ssshr(s)       ssmanage((s), SSSHRINK|SSSINGLE)
#define ssdelall()     ssmanage(NULL, SSDELETE)
#define ssshrall()     ssmanage(NULL, SSHRINK)
#define ssfgets(s, pn, fp)    ssfgetx(&(s), (pn), '\n', (fp))
#define ssfgetall(s, pn, fp)  ssfgetx(&(s), (pn), EOF, (fp))

ZCSTRCLS void ssmanage(char *, unsigned);
ZCSTRCLS char *sscpycatx(char **, const char *, size_t, unsigned);
ZCSTRCLS char *ssfgetx(char **, size_t *, int, FILE *fp);


#ifndef SSMINSIZ /* to override the block size, define it before inclusion */
#define SSMINSIZ 256 /* change this value to 1 for debugging */
#endif
#ifndef SSHASHBITS
#define SSHASHBITS 8
#endif
#define SSHASHSIZ  (1<<SSHASHBITS)  
#define SSOVERALLOC 1
#define sscalcsize_(n) (((n)/SSMINSIZ + 1) * SSMINSIZ) /* size for n nonblank characters */
#ifdef ZCOM_ERROR
#define sserror_ zcom_fatal
#else
/* print an error message and quit */
void sserror_(char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  fprintf(stderr, "fatal error: ");
  fprintf(stderr, fmt, args);
  va_end(args);
  exit(1);
}
#endif

struct ssheader{
  size_t size;
  size_t hashval;
  struct ssheader *next;
} ssbase_[SSHASHSIZ] = {{0u, 0u, NULL}};

/* we use the string address instead of that of the pointer 
 * to struct ssheader to compute the Hash value,
 * because the former is more frequently used in e.g. looking-up
 * */
static size_t sshashval_(const char *p)
{
  size_t val = (size_t)p * 1664525u + 1013904223u;
  return (val >> (sizeof(size_t)*8-SSHASHBITS)) & ((1<<SSHASHBITS)-1);
}

/* 
 * return the *previous* header to the one that associates with s
 * first locate the list from the Hash value, then enumerate the linked list.
 * */
static struct ssheader *sslistfind_(char *s)
{
  struct ssheader *hp;

  if (s == NULL)
    return NULL;
  for (hp = ssbase_ + sshashval_(s); hp->next != ssbase_; hp = hp->next)
    if ((char *)(hp->next + 1) == s)
      return hp;
  return NULL;
}

/* 
 * simply add the entry h at the begining of the list 
 * we do not accept a precalculated hash value, 
 * since realloc might have changed it
 * */
static struct ssheader *sslistadd_(struct ssheader *h)
{
  struct ssheader *head;

  head = ssbase_ + sshashval_( (char *)(h+1) );
  if (head->next == NULL) /* initialize the base */
    head->next = head;
  h->next = head->next;
  head->next = h;
  return head;
}

/* remove hp->next */
static void sslistremove_(struct ssheader *hp, int f)
{
  struct ssheader *h = hp->next;
  
  hp->next = h->next;
  if (f) free(h);
}

/* (re)allocate memory for (*php)->next, update list, return the new string
 * n is the number of nonempty characters, obtained e.g. from strlen().
 * create a new header if *php is NULL, in this case, the first character
 * of the string is '\0'
 * */
static char *ssresize_(struct ssheader **php, size_t n, unsigned flags)
{
  struct ssheader *h=NULL, *hp;
  size_t size;

  if (php == NULL)
    sserror_("NULL pointer to resize");
  
  /* we use the following if to assign hp and h, so the order is crucial */
  if ((hp=*php) == NULL || (h = hp->next)->size < n + 1 || !(flags & SSOVERALLOC)) {
    size = sscalcsize_(n);
    if (h == NULL || size != h->size) {
      /* since realloc will change the hash value of h
       * we have to remove the old entry first without free() 
       * hp->next will be freed by realloc */
      if (hp != NULL)
        sslistremove_(hp, 0);
      if ((h = realloc(h, sizeof(*h)+size)) == NULL) {
        sserror_("no memory for resizing\n");
        return NULL;
      }
      if (hp == NULL) /* clear the first byte if we start from nothing */
        *(char *)(h + 1) = '\0';  /* h + 1 is the beginning of the string */
      *php = hp = sslistadd_(h);
      hp->next->size = size;
    }
  }
  return (char *)(hp->next + 1);
}

static void ssmanage_low_(struct ssheader *hp, unsigned opt)
{
  if (opt == SSDELETE)
    sslistremove_(hp, 1);
  else if (opt == SSSHRINK)
    ssresize_(&hp, strlen((char *)(hp->next+1)), 0);
  else
    sserror_("unknown manage option");
}

/* delete a string, shrink memory, etc ... */
void ssmanage(char *s, unsigned flags)
{
  struct ssheader *hp,*head;
  unsigned opt = flags & 0xFF;
  size_t i;

  if (flags & SSSINGLE) {
    if (s == NULL || (hp = sslistfind_(s)) == NULL)
      return;
    ssmanage_low_(hp, opt);
  } else {
    for (i=0; i<SSHASHSIZ; i++)
      for (hp = head = ssbase_+i; hp->next && hp->next != head; hp = hp->next)
        /* we must not operate on h itself, which renders the iterator h invalid */
        ssmanage_low_(hp, opt);
  }
}

/* 
 * copy/cat t to *ps
 *
 * If (flags & SSCAT) == 0:
 * copy t to *ps, if ps is not NULL, and return the result
 * if ps or *ps is NULL, we return a string created from t
 *   *ps is set to the same value if ps is not NULL
 * otherwise, we update the record that corresponds to *ps
 *
 * minsize: to request a minimal size for the resulting buffer
 *
 * If flags & SSCAT:
 * append t after *ps. Equivalent to cpy if ps or *ps is NULL.
 * */
char *sscpycatx(char **ps, const char *t, size_t minsize, unsigned flags)
{
  struct ssheader *hp=NULL;
  size_t size=0u, sizes=0u;
  char *s=NULL, *p;

  /* both ps and *ps can be NULL, in which cases we leave hp as NULL */
  if (ps != NULL && (s=*ps) != NULL && (hp = sslistfind_(s)) == NULL)
    return NULL;
  if (t != NULL) 
    while (t[size]) /* compute the length of t */
      size++;
  if (flags & SSCAT) {
    if (s != NULL)  /* s is also NULL, if ps itself is NULL */
      while (s[sizes]) /* compute the length of s */
        sizes++;
    size += sizes;
  }  /* sizes is always 0 in case of copying */
  if (size < minsize)
    size = minsize;
  if ((s = ssresize_(&hp, size, SSOVERALLOC)) == NULL) /* change size */
    return NULL;
  if (t != NULL)
    for (p = s + sizes; (*p++ = *t++); ) /* copy/cat the string */
      ;
  if (ps != NULL)
    *ps = s;
  return s;
}

/* get a string *ps from file fp
 * *ps can be NULL, in which case memory is allocated
 * *pn is number of characters read (including '\n', but not the terminal null)
 * delim is the '\n' for reading a singe line
 * */
char *ssfgetx(char **ps, size_t *pn, int delim, FILE *fp)
{
  size_t n, max;
  int c;
  char *s;
  struct ssheader *hp;

  if (ps == NULL || fp == NULL)
    return NULL;
  if ((s=*ps) == NULL) /* allocate an initial buffer if *ps is NULL */
    if ((s = sscpycatx(ps, NULL, 0, 0u)) == NULL)
      return NULL;
  if ((hp = sslistfind_(s)) == NULL)
    return NULL;
  max = hp->next->size-1;
  for (n = 0; (c = fgetc(fp)) != EOF; ) {
    if (n+1 > max) { /* request space for n+1 nonblank characters */
      if ((*ps = s = ssresize_(&hp, n+1, SSOVERALLOC)) == NULL)
        return NULL;
      max = hp->next->size - 1;
    }
    s[n++] = (char)(unsigned char)c;
    if (c == delim) 
      break;
  }
  s[n] = '\0';
  if (pn != NULL)
    *pn = n;
  return (n > 0) ? s : NULL;
}
#endif /* ZCOM_SS__ */
#endif /* ZCOM_SS */


#ifdef  ZCOM_RNG
/* we have to define another macro ZCOM_RNG__ in order to avoid multiple inclusion
 * a single ZCOM_C__ won't do, because different module-set may be selected */
#ifndef ZCOM_RNG__
#define ZCOM_RNG__

#include <stdio.h>
#include <math.h>
#include <string.h>

/* a random number in [0,1) */
#define rnd0()        ((1.0/4294967296.0)*mtrand(0,0,NULL))  /* divide by 2^32 */
#define mtsave(fn)    mtrand(1,0,fn)
#define mtsetfile(fn) mtrand(2,0,fn)
#define mtfinish(fn)  mtrand(3,0,fn)

/* The following random number generator algorithm Mersenne Twister
   was developped by Makoto Matsumoto and Takuji Nishimura */

/* Period parameters */
#define N_MT 624
#define M_MT 397
#define UMASK_MT 0x80000000UL /* most significant w-r bits */
#define LMASK_MT 0x7fffffffUL /* least significant r bits */

/* The first argument defines the two possible actions
 * 0: to generate random number unsigned long
 *    for the first call, we try to load the state from file `fname',
 *    if it is unsuccessful, seed0, if unzero, is used to initialize the RNG
 * 1: save the current state to `fname'
 * 2: change the default file name
 * 3: save state and finish 
 * */
ZCSTRCLS unsigned long mtrand(int action, unsigned long seed0, const char *fname)
{
  static int mtindex_=-1; /* means it is not initialized */
  static unsigned long mt_[N_MT]; /* the array for the state vector  */
  static char *mtfname=NULL;

  unsigned long y;
  static unsigned long mag01[2]={0, 0x9908b0dfUL}; /* MATRIX_A */
  int k,nz=0;

  if (fname != NULL) /* if fname is given (not NULL), always use it */
    sscpy(mtfname, fname);
  else if (mtfname == NULL) /* otherwise if mtfname is missing, use the default */
    sscpy(mtfname, "MTSEED");

  if (action == 1 || action == 3) { /* to save the current state or finish */
    FILE *fp;
    if(mtindex_<0) return 1; /* RNG was never used, so it cannot be saved */
    if((fp=fopen(mtfname, "w"))==NULL){
      fprintf(stderr, "mtrand: cannot save state to %s.\n", mtfname);
      return 1;
    }
    fprintf(fp, "MTSEED\n%d\n", mtindex_);
    for(k=0; k<N_MT; k++)
      fprintf(fp, "%lu\n", mt_[k]);
    fclose(fp);
    if (action == 3) {
      ssdel(mtfname);
      mtfname = NULL;
    }
    return 0;
  }else if(action == 2){
    return 0;
  }else if(action != 0){
    fprintf(stderr, "mtrand: unknown action %d\n", action);
    return 0UL;
  }

  if(mtindex_ < 0){   /* initialize from the seed */
    FILE *fp;
    char s[64]="",err=1;

    /* here we try to initialize the array from file */
    if((fp=fopen(mtfname, "r")) != NULL){
      if(fgets(s, sizeof s, fp) == NULL){
        fprintf(stderr, "mtrand: cannot read the first line of %s, probably an empty file.\n", mtfname);
        goto CLOSEFILE;
      }

      if(strncmp(s, "MTSEED", 6) != 0){ /* to check the first line */
        fprintf(stderr, "mtrand: corrupted seed file.\n");
      }else{
        if( fscanf(fp, "%d", &mtindex_) != 1){
          fprintf(stderr, "mtrand: no index found.\n");
          goto CLOSEFILE;
        }
        if(mtindex_ < 0) mtindex_ = N_MT; /* request updating */
        for(k=0; k<N_MT; k++){
          if(EOF == fscanf(fp, "%lu", &mt_[k])) break;
          if(mt_[k]!=0) nz=1;
        }
        if(k != N_MT)
          fprintf(stderr, "mtrand: seed file is incomplete, got %d numbers.\n", k);
        else
          err=0; /* clear error mark */
      }
CLOSEFILE:
      fclose(fp);
    }
    if(!nz) err=1;
    if(err){
      unsigned long mtseed0=5489UL;
      if(seed0 != 0) mtseed0=seed0;
      mt_[0]= mtseed0 & 0xffffffffUL;
      for(k=1; k<N_MT; k++) {
        mt_[k] = (1812433253UL*(mt_[k-1]^(mt_[k-1]>>30))+k) & 0xffffffffUL;
      } /* the masking step is for 64-bit machines */
      mtindex_ = N_MT; /* request updating */
    }

  }

  if(mtindex_ >= N_MT) { /* generate N_MT words at one time */
    for(k=0; k<N_MT-M_MT; k++){
      y = (mt_[k]&UMASK_MT) | (mt_[k+1]&LMASK_MT);
      mt_[k] = mt_[k+M_MT] ^ (y>>1) ^ mag01[y&1UL];
    }
    for(; k<N_MT-1; k++){
      y = (mt_[k]&UMASK_MT) | (mt_[k+1]&LMASK_MT);
      mt_[k] = mt_[k+(M_MT-N_MT)] ^ (y>>1) ^ mag01[y&1UL];
    }
    y = (mt_[N_MT-1]&UMASK_MT) | (mt_[0]&LMASK_MT);
    mt_[N_MT-1] = mt_[M_MT-1] ^ (y>>1) ^ mag01[y&1UL];
    mtindex_ = 0;
  }

  y = mt_[mtindex_++];

  /* Tempering */
  y ^= (y>>11);
  y ^= (y<<7) & 0x9d2c5680UL;
  y ^= (y<<15) & 0xefc60000UL;
  y ^= (y>>18);

  return y;

#undef MTSETFNAME
}

#undef N_MT
#undef M_MT
#undef UMASK_MT
#undef LMASK_MT


/* Gaussian distribution with zero mean and unit variance,
   Algorithm adapted from Numerical recipes. */
ZCSTRCLS double grand0(void){
  const double tol=1e-15; /* 1e-13 allows fac of 1.0/(6e-8) */
  static int empty=1;
  static double save=0.0;
  double fac,r2,x,y;
  if(!empty) { empty=1; return save; }

  for(;;) {
    x=2.0*rnd0()-1.0;
    y=2.0*rnd0()-1.0;
    r2=x*x+y*y;
    /* to make sure it's inside a unit circle, and
       r2 can be a denominator */
    if(r2<1.0 && r2>tol) break;
  }
  fac=sqrt(-2.0*log(r2)/r2);
  empty=0;
  save=y*fac;
  return x*fac;
}

#endif /* ZCOM_RNG__ */
#endif /* ZCOM_RNG */

#ifdef  ZCOM_CFG
#ifndef ZCOM_CFG__
#define ZCOM_CFG__

/*
 * =======================================================================
 *
 * Configuration file
 *
 * ========================================================================
 *
 now cfgload should be replaced by cfgopen (more like fopen)

   USAGE:
  to load parameters use
  cfgget(fp, &var, "var_name", scanf_fmt);
  cfgget(fp, &arr_size, "arr_size", "%d");

  to save parameters, use fprintf :-)

for MPI, please only do it on the master node, then use MPI_Bcast
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct tagcfgdata_t{
  int n; /* number of lines */
  char**key,**value; /* key[i] is the key, value[i] is the value */
  char *buf; /* the whole configuration file */
  int dyn_alloc; /* whether the struct itself is dynamically allocated */
}cfgdata_t;

#define isspace_(c) isspace((unsigned char)(c))

/* load the whole configuration file into memory (return value); parse it to entries
   note: the order of cfg and filenm are exchanged */
ZCSTRCLS int cfgload(cfgdata_t *cfg, char *filenm)
{
  FILE *fp;
  long i,j;
  size_t size=0;
  char *p,*q,*lin;

  if((fp=fopen(filenm, "rb")) == NULL){
    fprintf(stderr,"cannot open the configuration file [%s]\n", filenm);
    return 1;
  }
  
  if (ssfgetall(cfg->buf, &size, fp) == NULL) {
    fprintf(stderr, "error reading file %s\n", filenm);
    return 4;
  }
  sscat(cfg->buf, "\n"); /* in case the file is not ended by a new line, we add one */
  fclose(fp);

#ifdef CFGDBG
  printf("the cfg is loaded, size=%d\n", size);
#endif

  /* count the number of lines for allocating the key-table */
  for (i = 0, cfg->n = 0; i<size; i++) {
    if (cfg->buf[i] == '\n' || cfg->buf[i] == '\r') {
      if (i>0 && cfg->buf[i-1]=='\\') {
        /* allow multiple-line splicing by replacing cr, lf with spaces
           parse should be aware of these additional spaces */
        cfg->buf[i-1]=' ';
        cfg->buf[i]=' ';
      } else {
        cfg->buf[i]='\0';
        (cfg->n)++;
      }

      for (j = i+1; j < size; j++) {
        /* we replace immediately followed cr & lf by spaces for
           efficiency (to avoid a large key table for blank lines) */
        if ( isspace_(cfg->buf[j]) ) {
          cfg->buf[j]=' ';
        } else {
          break;
        }
        /* note: parser should be insensitive to leading spaces */
      }
#ifdef CFGDBG
      if (j-1 >= i+1) printf("j=%d to %d are replaced by spaces\n", i+1,j-1);
#endif
    }
  }
#ifdef CFGDBG
  printf("# of lines: %d\n", cfg->n);
#endif

  cfg->key   = calloc(cfg->n, sizeof(char*));
  cfg->value = calloc(cfg->n, sizeof(char*));

  /* load lines into the keytable, not parsed yet */
  for(p=q=cfg->buf,j=0,i=0; i < size; i++){
    if(cfg->buf[i]=='\0'){
      cfg->key[j] = p;
      j++;
      p=cfg->buf+i+1;
      /* we may still have spaces left over, but no need to continue */
      if(j>cfg->n) break;
    }
  }
  cfg->n=j;
#ifdef CFGDBG
  fprintf(stderr, "load %d lines.\n", cfg->n);
#endif

  /* now parse lines: separate values from keys */
  for(j=0; j<cfg->n; j++){
    lin=cfg->key[j];

    /* remove the the leading spaces */
    for(; *lin && isspace_(*lin); lin++) ;
    cfg->key[j]=lin;
    /* skip a blank or comment line */
    if(lin[0]=='\0' || strchr("#%!;", lin[0])!=NULL){
      cfg->key[j]=NULL;
      continue;
    }

    /* remove trailing space and ';' */
    for(q =lin+strlen(lin)-1;
        q>=lin && (isspace_(*q)||*q==';'); q--) *q='\0';

    if((q=strchr(lin, '=')) == NULL){ /* skip a line without '=' */
      cfg->key[j]=NULL;
      continue;
    }

    /* find the end of key --> 'q' */
    *q='\0';
    p=q+1;
    for(--q; isspace_(*q); q--) *q='\0';
    for(; (*p) && isspace_(*p); p++) ; /* skip leading space, 'p' -> value */
    cfg->value[j]=p;
  }

#ifdef CFGDBG
  for(j=0; j<cfg->n; j++){ if(cfg->key[j]) printf("key=%s, value=%s\n",cfg->key[j], cfg->value[j]); }
  printf("%d lines\n",cfg->n); getchar();
#endif

  return 0;
}

#undef isspace_

/* a wrapper of cfgload to make it more like fopen */
ZCSTRCLS cfgdata_t *cfgopen(char *filenm){
  cfgdata_t *cfg;
  if((cfg=calloc(1,sizeof(cfgdata_t)))==NULL){
    fprintf(stderr, "cannot allocate space for cfgdata_t.\n");
    return NULL;
  }
  if(cfgload(cfg,filenm) != 0){
    free(cfg);
    return NULL;
  }
  cfg->dyn_alloc=1;
  return cfg;
}


ZCSTRCLS void cfgclose(cfgdata_t *cfg){
  free(cfg->value); cfg->value=NULL;
  free(cfg->key);   cfg->key=NULL;
  ssdel(cfg->buf);  cfg->buf=NULL;
  if(cfg->dyn_alloc) free(cfg);
}


/* read the value of a given variable from the current configuration file,
   the name of variable is given by `key',
   if the key is matched, its value is saved to `*var' through sscanf.
   if the function succeeds, it returns 0.
   A comment line in the configuration file starts with '#', '%' or '!'.
   In case fmt is "%s", a duplicate of the string will be assigned to (*var)
   */
ZCSTRCLS int cfgget(cfgdata_t *cfg, void *var, const char *key, const char* fmt){
  int j;

  if (cfg->key == NULL || var==NULL || key==NULL || fmt==NULL) {
    fprintf(stderr, "cfgget: NULL pointer.\n");
    return 1;
  }

  for (j = 0; j < cfg->n; j++) 
    if (cfg->key[j] != NULL && strcmp(cfg->key[j], key) == 0) { 
      if (fmt[0]=='%' && fmt[1]=='s') { /* string case */
        *(char **)var = ssnew(cfg->value[j]); /* make a copy and return */
        return 0;
      } else { /* use sscanf for other cases, like int,float,... */
        if (EOF == sscanf(cfg->value[j], fmt, var))
          return 2; /* input error */
        else
          return 0;
      }
    }
  return 1; /* no match */
}

#ifdef ZCOM_CFG_LEGACY
/* read the value of a given array variable from the current configuration file,
   the size of each array item is given by `itemsize',
   the name of the array starts with `key' followed by an index
   e.g. if key="arr", entries in configuration files are `arr0', `arr1', ...
   indices are from i0 to iN-1;
   if the function succeeds, it returns 0. */
ZCSTRCLS int cfggetarr(cfgdata_t *cfg, void const *varr, size_t itemsize, const char *key, const char* fmt, int i0, int iN){
  int i;
  char *var, itemname[128];
  var=(char*)varr;
  if(strlen(key)>sizeof(itemname)-16){
    fprintf(stderr, "key name is too long\n");
    return 1;
  }
  for(i=i0; i<iN; i++){
    sprintf(itemname, "%s%d", key, i);
    cfgget(cfg, var, itemname, fmt);
    var+=itemsize;
  }
  return 0;
}
#endif

#ifdef ZCOM_CFG_TST
int main(void){
  int nr, cnt, i;
  float tmin, tmax;
  float arr[10];
  char*p=NULL;
  cfgdata_t *cfg;
  if((cfg=cfgopen("test.cfg")) == NULL){
    printf("error reading\n");
    return 1;
  }
  cfgget(cfg, &nr, "nrtemp", "%d");
  cfgget(cfg, &tmin, "tmin", "%f");
  cfgget(cfg, &tmax, "tmax", "%f");
  cfgget(cfg, &cnt, "arrcnt", "%d");
  cfggetarr(cfg, arr, sizeof(float), "arr", "%f", 0, cnt);
  cfgget(cfg, &p, "scode", "%s");
  cfgclose(cfg);

  printf("nr=%d, (%g,%g), cnt=%d\n", nr,tmin,tmax,cnt);
  for(i=0; i<cnt; i++){
    printf("arr[%d]=%g\n", i, arr[i]);
  }
  printf("scode=\"%s\"\n", p);

  return 0;
}
#endif

#endif /* ZCOM_CFG__ */
#endif /* ZCOM_CFG */

#ifdef  ZCOM_TRACE
#ifndef ZCOM_TRACE__
#define ZCOM_TRACE__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#ifdef ZCHAVEVAM
/* we use variadic macros if possible */
#define wtrace(fmt, ...)     wtrace_x(0, fmt, ## __VA_ARGS__)
#define wtrace_buf(fmt, ...) wtrace_x(1, fmt, ## __VA_ARGS__)
/* Buffered (flags=0) or unbuffered (flags=1) output
 *
 * pass NULL to `fmt' to finish up
 * In addition, we support simple commands in `fmt'
 * If `fmt' starts with %@, we start the command mode
 * if fmt="%@filename=%s", the next argument is "tr.dat",
 *    then the trace file becomes tr.dat
 * if fmt="%@freq=%d", the next argument is 100,
 *    then we refresh every 100 calls.
 * */
ZCSTRCLS int wtrace_x(int wtrace_flags_, const char *fmt, ...)

#else

/* if variadic macros are unavailable, we only define the buffered version */
#define wtrace wtrace_buf
static int wtrace_flags_=1;
ZCSTRCLS int wtrace_buf(const char *fmt, ...)
#endif
{
  const int maxlen=1024;
  static int verbose=1;
  static int flush_interval=1000;
  static int cnt=0, once=0, i;
  static char *buf,*msg=NULL;
  static FILE *fp=NULL;
  static char *fname=NULL, mode[8]="w";
  va_list args;

  if (fname == NULL) /* set the default file name */ 
    fname = ssnew("TRACE");

  /* to finish up */
  if(fmt == NULL) goto NORMAL;

  /* start the command mode if the format string start with "%@"
   * the command mode allows setting parameters */
  if (fmt[0] == '%' && fmt[1] == '@') {
    const char *cmd=fmt+2, *p;
    /* we no longer allow commands after the tracing starts */
    if(once) goto NORMAL;

    if((p=strchr(fmt, '='))==NULL) goto NORMAL;
    va_start(args, fmt);
    if(strncmp(cmd, "filename", 8)==0){
      p=va_arg(args, const char *);
      if (p != NULL) {
        sscpy(fname, p);
        if (verbose)
          fprintf(stderr, "The trace file is now %s\n", fname);
      }
    }else if(strncmp(cmd, "freq", 4) == 0){
      flush_interval=va_arg(args, int);
    }else if(strncmp(cmd, "verbose", 7) == 0){
      verbose=va_arg(args,int);
    }else{ /* unknown command */
      goto NORMAL;
    }
    va_end(args);
    return 0;
  }

NORMAL:
  if(wtrace_flags_ == 0){ /* unbuffered version */
    if(cnt==0){
      if(once) mode[0]='a'; else mode[0]='w';
      if( (fp=fopen(fname, mode)) == NULL ){
        fprintf(stderr, "cannot write file %s with mode %s\n", fname, mode);
        return 1;
      }
    }

    if(fmt != NULL){
      va_start(args, fmt);
      vfprintf(fp, fmt, args);
      va_end(args);
    }

    if((++cnt == 1000) || fmt==NULL){
      cnt=0;
      fclose(fp);
    }
    once=1;
    return 0;
  }

  if(cnt==0){
    /* we allocate msg and buf together, because msg may be (though unlikely)
     * involved in an unsafe vsprintf(), which might cause overflow.
     * in this case, the continuous space can help prevent a system crash,
     * although buf can be corrupted */
    msg=calloc(maxlen*(flush_interval+1), 1);
    if(msg==NULL){
      fprintf(stderr, "cannot allocate buffer.\n");
      return 1;
    }
    msg[0]='\0';
    buf=msg+maxlen;
    buf[0]='\0';
  }

  if(msg==NULL) return 1;

  if(fmt != NULL){
    if(strlen(fmt) >= (size_t)maxlen){
      fprintf(stderr, "the format string is too long.\n");
      return 1;
    }

    va_start(args, fmt);
#ifdef _MSC_VER
  #ifndef vsnprintf
    #define vsnprintf _vsnprintf
  #endif
#endif

#if (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__xlC__) || defined(_MSC_VER) )
    i=vsnprintf(msg, maxlen, fmt, args);
#else
    /* we use vsprintf if vsnprintf is unavailable */
    i=vsprintf(msg, fmt, args);
#endif
    va_end(args);
    if(i>=maxlen){
      fprintf(stderr, "the message is too long.\n");
      return 1;
    }
    strcat(buf, msg);
  }

  if((++cnt % flush_interval == 0) || fmt==NULL){
    mode[0] = (char)(once ? 'a' : 'w');
    if ( (fp=fopen(fname, mode)) == NULL ) {
      fprintf(stderr, "cannot write file %s with mode %s\n", fname, mode);
      return 1;
    }
    fputs(buf, fp);
    buf[0]='\0';
    fclose(fp);

    if (fmt == NULL) { /* finishing up */
      if (msg != NULL) { 
        free(msg); 
        msg = NULL; 
      }
      if (fname != NULL) {
        ssdel(fname);
        fname = NULL;
      }
      cnt = 0;
      once = 0;
    }else{  /* subsequent calls will append instead of write the trace file */
      once=1;
    }
  }
  return 0;
}
#endif /* ZCOM_TRACE__ */
#endif /* ZCOM_TRACE */


#ifdef  ZCOM_DIST
#ifndef ZCOM_DIST__
#define ZCOM_DIST__

/* to normalize and write histograms to file
 pointer 'h' (to be converted from 1D or 2D array) contains 'rows' histograms,
 each contains 'cols' entries, from 'base' to 'base+inc*cols' */
#define WD_ADDAHALF   0x0001
#define WD_KEEPEDGE   0x0002
#define WD_NOZEROES   0x0004
#define wdist(h,rows,cols,base,inc,fname)  wdistex(h,rows,cols,base,inc,WD_ADDAHALF,fname)

ZCSTRCLS int wdistex(double *h, int rows, int cols, double base, double inc, int flag, char *fname)
{
  char *filename;
  FILE *fp;
  int i,j,imax,imin;
  double sum,*p,delta;

  filename = ssnew((fname != NULL) ? fname : "HIST");
  
  if ((fp=fopen(filename, "w")) == NULL) {
    printf("cannot write history file [%s].\n", filename);
    return 1;
  }
  delta = (flag & WD_ADDAHALF) ? 0.5 : 0;

  for(j=0; j<rows; j++){
    p=h+j*cols;

    if(flag&WD_KEEPEDGE){
      imin=0;
      imax=cols;
    }else{
      for(i=cols-1; i>=0; i--) if(p[i]>0) break;
      imax=i+1;
      if(imax == 0) continue;

      for(i=0; i<imax; i++) if(p[i]>0) break;
      imin=i;
    }

    for(sum=0,i=imin; i<imax; i++) sum += p[i];
    sum *= inc;
    if(fabs(sum)<1e-6) sum=1;

    for(i=imin; i<imax; i++){
      if((flag&WD_NOZEROES) && p[i]<1e-6) continue;
      fprintf(fp,"%g %20.14E %d\n", base+(i+delta)*inc, p[i]/sum, j);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  
  ssdel(filename);
  return 0;
}

#endif /* ZCOM_DIST__ */
#endif /* ZCOM_DIST */


#ifdef  ZCOM_LOG
#ifndef ZCOM_LOG__
#define ZCOM_LOG__
/*
 * =======================================================================
 *
 * LOG file routines
 *
 * ========================================================================
 */
typedef struct tag_logfile_t{
  FILE *fp;
  char *fname;
  int flag;
}logfile_t;

#define LOG_WRITESCREEN  0x01
#define LOG_FLUSHAFTER   0x02
#define LOG_NOWRITEFILE  0x10
ZCSTRCLS logfile_t* logopen(char *filenm)
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
  log->fname = ssnew(filenm);
  log->flag=0;
  return log;
}

ZCSTRCLS int logprintf(logfile_t *log, char *fmt, ...)
{
  va_list args;
  if(log==NULL) return 1;

  if(log->fp==NULL) log->fp = fopen(log->fname, "w");
  if(log->fp==NULL){
    fprintf(stderr, "log [%s] cannot be opened.\n", log->fname);
    return 1;
  }
  if((log->flag&LOG_NOWRITEFILE) == 0){
    va_start(args, fmt);
    vfprintf(log->fp, fmt, args);
    va_end(args);
  }
  if(log->flag&LOG_WRITESCREEN){
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
  if(log->flag&LOG_FLUSHAFTER) fflush(log->fp);
  return 0;
}

/* close & reopen log file to make sure that stuff is written to disk */
ZCSTRCLS int loghardflush(logfile_t *log)
{
  if(log->fp==NULL || log->fname==NULL) return 1;
  fclose(log->fp);
  if((log->fp=fopen(log->fname, "a"))==NULL){
    fprintf(stderr, "cannot reopen the log file [%s].\n",log->fname);
    return 1;
  }
  return 0;
}

ZCSTRCLS void logclose(logfile_t *log) 
{
  if (log == NULL) 
    return;
  if (log->fp != NULL) {
    fclose(log->fp);
    log->fp = NULL;
  }
  if (log->fname != NULL) {
    ssdel(log->fname);
    log->fname = NULL;
  }
  free(log);
}

#ifdef LOGTST
int main(void)
{
  logfile_t *mylog;
  char msg[256];
  if((mylog=logopen("my.log")) == NULL) return 1;
  printf("please write something: ");
  logprintf(mylog, "the input is [%s]\n", fgets(msg, sizeof msg, stdin));
  loghardflush(mylog);
  printf("log file is hard flushed. please check...");
  getchar();
  logprintf(mylog, "finished.\n");
  logclose(mylog);
  return 0;
}
#endif

#endif /* ZCOM_LOG__ */
#endif /* ZCOM_LOG */



#ifdef  ZCOM_LU
#ifndef ZCOM_LU__
#define ZCOM_LU__
/* LU decomposition part : */

/* solve A x = b by L U decomposition */
ZCSTRCLS int lusolve(double *a, double *b, int n){
  int i,j,k,imax=0;
  double x,max,sum;
  const double mintol=1e-16; /* absolute minimal value for a pivot */

  for (i=0;i<n;i++) {  // normalize each row
    for(max=0.0,j=0;j<n;j++)
      if ((x=fabs(a[i*n+j])) > max) max=x;
    if(max<mintol){
      printf("lusolve: coefficients too small.\n");
      return 1;
    }
    for(x=1.0/max,j=0; j<n; j++) a[i*n+j]*=x;
    b[i]*=x;
  }

  // step 1: A = L U, column by column
  for(j=0;j<n;j++) {
    // matrix U
    for(i=0;i<j;i++){
      sum=a[i*n+j];
      for(k=0;k<i;k++) sum-=a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
    }

    // matrix L, diagonal of L are 1
    max=0.0;
    for (i=j;i<n;i++) {
      sum=a[i*n+j];
      for(k=0;k<j;k++) sum-=a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
      if((x=fabs(sum)) >= max) { max=x; imax=i; }
    }

    if (j != imax) { // swap the pivot row with the jth row
      for (k=0;k<n;k++) { x=a[imax*n+k]; a[imax*n+k]=a[j*n+k]; a[j*n+k]=x; }
      x=b[imax]; b[imax]=b[j]; b[j]=x;
    }
    if(fabs(a[j*n+j])<mintol){
      printf("lusolve: cannot pick a pivot for %d for %dD.\n", j, n);
      return 2;
    }
    // divide by the pivot element, for the L matrix
    if(j!=n-1) for(x=1.0/a[j*n+j],i=j+1;i<n;i++) a[i*n+j]*=x;
  }

  // step2: solve the equation L U x = b
  for(i=0;i<n;i++) { // L y = b
    x=b[i];
    for(j=0;j<i;j++) x-=a[i*n+j]*b[j];
    b[i]=x;
  }
  for(i=n-1;i>=0;i--) { // U x = y.
    x=b[i];
    for(j=i+1;j<n;j++) x -= a[i*n+j]*b[j];
    b[i]=x/a[i*n+i];
  }
  return 0;
}


/*
  LU decomposition of an n x n matrix 'a',
  The diagonal belongs to U, the diagonal of L is 1's
  idx[] as output can have duplicated indices, but this is readable by
*/
ZCSTRCLS int ludcmp(double *a, int *idx, int n){
  int i,j,k,imax;
  double x,max,sum, *scal;
  const double mintol=1e-16; /* absolute minimal value for a pivot */

  if((scal=(double *)malloc(sizeof(double)*n)) == NULL){
    printf("cannot allocate vector\n");
    return -1;
  }

  for(i=0;i<n;i++) {  // normalize each row
    for(max=0,j=0;j<n;j++)
      if ((x=fabs(a[i*n+j])) > max) max=x;
    if(max<mintol){
      printf("ludcmp: zero equation.\n");
      free(scal);
      return 1;
    }
    scal[i]=1.0/max;
    //printf("scal[%d]=%g\n", i, scal[i]);
  }

  // decompose  A = L U, column by column, first U, next L
  for(j=0;j<n;j++) {
    // matrix U
    for(i=0;i<j;i++){
      sum=a[i*n+j];
      for(k=0;k<i;k++) sum-=a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
    }

    // matrix L, diagonal of L are 1
    max=0.0; imax=j;
    for(i=j; i<n; i++) {
      sum=a[i*n+j];
      for(k=0; k<j; k++) sum-=a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
      // choosing the biggest element as the pivot
      if((x=scal[i]*fabs(sum)) >= 1.001*max) { max=x; imax=i; }
    }
    if (j != imax) { // swap the pivot row and the jth row
      for(k=0; k<n; k++) { x=a[imax*n+k]; a[imax*n+k]=a[j*n+k]; a[j*n+k]=x; }
      scal[imax]=scal[j]; // don't care scal[j], no use in the future
    }
    idx[j]=imax;
    //printf("j=%d, imax=%d, max=%g\n", j, imax, max);
    if(fabs(a[j*n+j])<mintol){
      printf("ludcmp: cannot pick a pivot for %d.\n", j);
      free(scal);
      return 2;
    }
    // divide by the pivot element, for the L matrix
    if(j!=n-1) for(x=1.0/a[j*n+j],i=j+1;i<n;i++) a[i*n+j]*=x;
  }
  free(scal);
  return 0;
}

/* backward substitution for L U x = b */
ZCSTRCLS int lubksb(double *a, double *b, int *idx, int n){
  int i, j, id;
  double x;

  for(i=0; i<n; i++) // unscramble indices
    if((id=idx[i]) != i) { x=b[id]; b[id]=b[i]; b[i]=x; }

  for(i=0; i<n; i++) { // L y = b
    for(x=b[i],j=0; j<i; j++) x-=a[i*n+j]*b[j];
    b[i]=x;
  }
  for(i=n-1; i>=0; i--) { // U x = y.
    x=b[i];
    for(j=i+1; j<n; j++) x-=a[i*n+j]*b[j];
    b[i]=x/a[i*n+i];
  }
  return 0;
}
#ifdef LUTST
#define N 3
int main(void){
  double A[N][N]={{1,0.5,1.0/3},{0.5,1.0/3,0.25},{1.0/3,.25,.2}}, a[N][N];
  double b[N]={1,2,3}, x[N];
  int i,j,index[N];

  memcpy(a, A, sizeof(double)*N*N);
  memcpy(x, b, sizeof(double)*N);
  lusolve((double *)a, x, N);
  for(i=0; i<N; i++) printf("%d: %g\n", i, x[i]);

  memcpy(a, A, sizeof(double)*N*N);
  memcpy(x, b, sizeof(double)*N);
  ludcmp((double *)a, index, N);
  for(i=0; i<N; i++){
    printf("%d: ", index[i]);
    for(j=0; j<N; j++){
      printf("%10.6f ", a[i][j]);
    }
    printf("\n");
  }

  lubksb(a, x, index, N);
  for(i=0; i<N; i++) printf("%d: %g\n", i, x[i]);

  return 0;
}
#endif /* LUTST */

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
  size_t i,j;
  int cs, ct, docase,doupper;

  docase=(flags&ZCOM_DOCASE); /* do case conversion */
  doupper=(flags&ZCOM_UPPER);

  if(flags&ZCOM_STRCMP){ /* comparison, size_s is ignored */
    if(s==NULL||t==NULL) return 0;
    for(i=0; ; i++){
      cs=s[i];
      ct=t[i];
      if(docase){
        cs=tolower((unsigned char)cs);
        ct=tolower((unsigned char)ct);
      }
      if(cs==0 || ct==0 || cs != ct) break;
    }
    return cs-ct;
  }else if(flags & (ZCOM_STRCPY|ZCOM_STRCAT)){ /* copying and */
    if(size_s==0 || s==NULL || t==NULL) return 1;
    /* t[size_s-1] should be '\0' */
    i=0;
    if(flags&ZCOM_STRCAT) while(s[i]) i++;
    for(j=0; i<size_s-1; i++,j++){
      ct=t[j];
      if(docase && (ct!=0)){
        ct=(unsigned char)(char)ct;
        if(doupper){
          cs=toupper(ct);
        }else{
          cs=tolower(ct);
        }
      }else{
        cs=ct;
      }
      s[i]=(char)cs;
      if(ct==0) break;
    }
    if(i==size_s-1) s[i]='\0';
    return (t[j]!='\0');
  }else{ /* unknown flag */
    fprintf(stderr, "zcom_strconv: invalid flags=%#x.\n", flags);
    return 0;
  }
}

#endif /* ZCOM_STR__ */
#endif /* ZCOM_STR */


#ifdef  ZCOM_ZT
#ifndef ZCOM_ZT__
#define ZCOM_ZT__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#define rzt(lnz,beta,cnt,erg0,max,flag,fname) rztx(lnz,beta,NULL,NULL,cnt,erg0,max,flag,fname)

/* write ln(Z) versus temperature */
ZCSTRCLS int wztf(double *lnz, double *beta, double *hist, int rows, int cols,
         double base, double inc, double erg0, const char *fname){
  char *filename;
  FILE *fp;
  int i,j;
  double c,x,x2,*h;

  filename = ssnew((fname != NULL) ? fname : "ZT");
  if((fp=fopen(filename, "w"))==NULL || beta==NULL){
    fprintf(stderr, "cannot write file [%s].\n", fname);
    return 1;
  }
  for(j=0; j<rows; j++){
    x=x2=c=0;
    if(hist){
      h=hist+j*cols;

      for(i=0; i<cols; i++){
        double y,E;
        if((y=h[i]) > 1e-6){
          c+=y; y*=(E=base+i*inc);
          x+=y; y*=E;
          x2+=y;
        }
      }
      if(c > 1e-6){
        x/=c;
        x2/=c;
        x2-=x*x;
      }
    }
    fprintf(fp, "%g %16.12f %16.12f %16.12f %g\n", 1.0/beta[j],
      ((lnz!=NULL)?(lnz[j]-beta[j]*erg0):0), x, (beta[j]*beta[j])*x2, c);
  }
  fclose(fp);
  ssdel(filename);
  return 0;
}

ZCSTRCLS int rztx(double *lnz, double *beta, double *erg, double *cv, int *cnt, double erg0, int max, int flag, char *fname){
  double T;
  FILE *fp;
  char s[1024];
  double enrg,hcap;

  if((fp=fopen(fname, "r"))==NULL){
    printf("cannot read file [%s]", fname);
    return 1;
  }
  for((*cnt)=0; fgets(s, sizeof s, fp); ){
    if((*cnt)>=max) break;
    if(s[0] == '#') continue; // comment line
    sscanf(s, "%lf%lf%lf%lf\n", &T, &lnz[*cnt], &enrg, &hcap);
    if(erg != NULL) erg[*cnt]=enrg;
    if(cv != NULL) cv[*cnt]=hcap;
    if(flag==0){
      beta[*cnt]=1.0/T;
    }else if(flag==1){ // this option checks if the beta matches
      if(fabs(1.0/beta[*cnt]-T)>1e-3){
        fprintf(stderr, "file is not compatible with current setup.\n"
        "line %d: expect %g, got %g\n", *cnt, 1.0/beta[*cnt], T);
        return 1;
      }
    }else{
      fprintf(stderr, "unknown flag %d.\n", flag);
      return 2;
    }
    lnz[*cnt] += erg0/T;
    (*cnt)++;
  }
  fclose(fp);
  return 0;
}

int vfscanf(FILE *fp, const char *fmt, va_list ap);

/* read and write simple data to file */
ZCSTRCLS int rwonce(const char *flag, const char *fname, const char *fmt, ...){
  FILE *fp;
  va_list args;

  if((fp=fopen(fname, flag))==0){
    printf("cannot open file [%s]\n", fname);
    return 1;
  }

  if(flag[0]=='w'){
    va_start(args, fmt);
    vfprintf(fp, fmt, args);
    va_end(args);
  }else{

#ifdef _MSC_VER
    void *ptr;
    char *p,*q,*format,ch;
    int i;

    format = ssnew(fmt);

    /* to seek the first conversion */
    for(p=format; *p; p++) if(*p=='%'){ if(p[1]=='%' || p[1]=='*') p++; else break;}
    if(*p){
      va_start(args, fmt);
      ptr=va_arg(args, void *);
      do{
        for(q=p+1; *q; q++) if(*q=='%'){ if(q[1]=='%' || q[1]=='*') q++; else break;}
        if( (ch=*q) != 0) *q='\0';
        i=fscanf(fp, p, ptr);
        if(i && ch) ptr=va_arg(args, void *); /* next conversion */
        (*(p=q)) = ch;
      }while(*q);
      va_end(args);
    }else{
      fscanf(fp, format); /* no conversion */
    }
    ssdel(format);

#else
    va_start(args, fmt);
    if( vfscanf(fp, fmt, args) > 0)  /* empty */;
    va_end(args);
#endif

  }
  fclose(fp);
  return 0;
}

// return log(exp(a)+exp(b))
#define LOGBIG 40.0  // error is 4e-18
static __inline double logadd(double a, double b){
  double c;
  if(a>b) return (  ((c=a-b)>LOGBIG) ? a : (a+log(1+exp(-c)))  );
  else    return (  ((c=b-a)>LOGBIG) ? b : (b+log(1+exp(-c)))  );
}

#define RZ_ADDAHALF   0x0001
#define RZ_OVERWRITE  0x0002

// refine lnZ through iteration
ZCSTRCLS int refinelnZ(double beta[], double X[], double *lnZbeta0, double *h, double *lng,
  int Tcnt, int Ecnt, double Ebase, double Einc, int flag,
  char *filename, char *DOSname, double err, int repmax)
{
  FILE *fp;
  int i,j,icnt;
  double Ehalf,E,c,diff=0;
  double mZ0=0,num,den;
  static double *lnZ, *lnN, *lnm;
  int rep;

#define LOG0    -1e10 // exp(LOG0) is almost 0

  lnZ=malloc(sizeof(double)*Tcnt);
  lnN=malloc(sizeof(double)*Tcnt);
  if(lnZ==NULL || lnN==NULL){
    printf("cannot allocate space for lnZ or lnN.\n");
    return 1;
  }

  // shift the partition function
  for(j=1; j<Tcnt; j++) lnZ[j]=X[j]-X[0];
  lnZ[0]=0;

  // determine the energy range
  for(icnt=0, j=0; j<Tcnt; j++){
    for(i=icnt; i<Ecnt; i++)
      if(h[j*Ecnt+i]>0.1 && i>icnt) icnt=i;
  }
  icnt++;
  lnm=malloc(sizeof(double)*icnt);
  if(lnm==NULL){
    printf("cannot allocate space for lnm.\n");
    free(lnZ); free(lnN);
    return 1;
  }
  for(i=0; i<icnt; i++) lnm[i]=0;

  // logN[j] is the total visits to temperature j
  for(j=0; j<Tcnt; j++){
    for(c=i=0; i<icnt; i++){
      c+=num=h[j*Ecnt+i];
      lnm[i]+=num;
    }
    lnN[j]=log(c);
  }
  for(i=0; i<icnt; i++)
    lnm[i]=((lnm[i]>0.1)?log(lnm[i]):LOG0);

  if(flag&RZ_ADDAHALF) Ehalf=Einc*0.5; else Ehalf=0;

  // repeat until convergence
  for(rep=1; rep<repmax; rep++){
    // calculate the logarithm of the density of states
    for(i=0; i<icnt; i++){
      E=(Ebase+i*Einc+Ehalf);
      num=0;
      den=LOG0; // so that exp(den)=0;
      for(j=0; j<Tcnt; j++)
        den=logadd(den, lnN[j]-beta[j]*E-lnZ[j]);
      lng[i] = ((lnm[i]<LOG0+0.1)?LOG0:(lnm[i]-den));
    }

    diff=0;
    for(j=0; j<Tcnt; j++){
      // calculate the new partition function
      c=LOG0; // c=log(Z[t]);
      for(i=0; i<icnt; i++)
        c=logadd(c, lng[i]-beta[j]*(Ebase+i*Einc+Ehalf));

      if(j==0) mZ0=c;
      else{
        c-=mZ0;
        if(fabs(lnZ[j]-c)>diff) diff=fabs(lnZ[j]-c);
        lnZ[j]=c;
      }
    }
    if(diff<err) break;
  }
  if(diff>err) printf("Warning: failed to converge! diff=%g\n", diff);
  else printf("converged after %d iterations.\n", rep);

  // write the density of states
  if((fp=fopen(DOSname, "w"))==NULL){
    printf("cannot open DOS\n");
  }else{
    for(i=0; i<icnt; i++)
      if(lng[i] > LOG0)
        fprintf(fp, "%g %16.10f\n", Ebase+i*Einc+Ehalf, lng[i]);
    fclose(fp);
  }

  // write the refined partition function
  if((fp=fopen(filename, "w"))==NULL){
    printf("cannot open [%s]\n", filename);
  }else{
    for(j=0; j<Tcnt; j++)
      fprintf(fp, "%g %16.10f\n", 1.0/beta[j], lnZ[j]);
    fclose(fp);
  }

  if(lnZbeta0){
    for(c=LOG0,i=0; i<icnt; i++) c=logadd(c, lng[i]);
    *lnZbeta0=c;
  }

  if(flag&RZ_OVERWRITE){
    for(j=0; j<Tcnt; j++) X[j]=lnZ[j];
  }

  free(lnZ); free(lnN); free(lnm);
  return 0;
}

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

ZCINLINE void rv3_zero(real *x){x[0]=0.0f; x[1]=0.0f; x[2]=0.0f; }
ZCINLINE void rv3_copy(real *x, const real *src){x[0]=src[0]; x[1]=src[1]; x[2]=src[2]; }

ZCINLINE real rv3_sqr (const real *x){return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];}
ZCINLINE real rv3_norm(const real *x){return (real)sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);}

ZCINLINE real *rv3_normalize(real *x){
  real r=rv3_norm(x);
  if(r>0.0){
    r=1.0f/r;
    x[0]*=r;
    x[1]*=r;
    x[2]*=r;
  }
  return x;
}

/* if x == y, try to use sqr */
ZCINLINE real rv3_dot(const real *x, const real *y)
{
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

ZCINLINE real *rv3_cross(real *ZCRESTRICT z, const real *x, const real *y){
  z[0]=x[1]*y[2]-x[2]*y[1];
  z[1]=x[2]*y[0]-x[0]*y[2];
  z[2]=x[0]*y[1]-x[1]*y[0];
  return z;
}

ZCINLINE real *rv3_neg(real *x){
  x[0]-=x[0];
  x[1]-=x[1];
  x[2]-=x[2];
  return x;
}

ZCINLINE real *rv3_neg2(real *nx, const real *x){
  nx[0]=-x[0];
  nx[1]=-x[1];
  nx[2]=-x[2];
  return nx;
}

ZCINLINE real *rv3_inc(real *x, const real *dx){
  x[0]+=dx[0];
  x[1]+=dx[1];
  x[2]+=dx[2];
  return x;
}
ZCINLINE real *rv3_dec(real *x, const real *dx){
  x[0]-=dx[0];
  x[1]-=dx[1];
  x[2]-=dx[2];
  return x;
}
ZCINLINE real *rv3_sinc(real *x, real *dx, real s){
  x[0]+=s*dx[0];
  x[1]+=s*dx[1];
  x[2]+=s*dx[2];
  return x;
}
ZCINLINE real *rv3_smul(real *x, real s){
  x[0]*=s;
  x[1]*=s;
  x[2]*=s;
  return x;
}

/* if y == x, just use smul */
ZCINLINE real *rv3_smul2(real * ZCRESTRICT y, const real *x, real s){
  y[0]=x[0]*s;
  y[1]=x[1]*s;
  y[2]=x[2]*s;
  return y;
}

/* for in-place difference use rv3_dec */
ZCINLINE real *rv3_diff(real * ZCRESTRICT diff, const real *a, const real *b){
  diff[0]=a[0]-b[0];
  diff[1]=a[1]-b[1];
  diff[2]=a[2]-b[2];
  return diff;
}

/* sum = a+b, for in-place addition use rv3_inc */
ZCINLINE real *rv3_sum2(real * ZCRESTRICT sum, const real *a, const real *b){
  sum[0]=a[0]+b[0];
  sum[1]=a[1]+b[1];
  sum[2]=a[2]+b[2];
  return sum;
}

/* sum = -a-b */
ZCINLINE real *rv3_nsum2(real *sum, const real *a, const real *b){
  sum[0]=-a[0]-b[0];
  sum[1]=-a[1]-b[1];
  sum[2]=-a[2]-b[2];
  return sum;
}

ZCINLINE real *rv3_lincomb2(real *sum, const real *a, const real *b, real s1, real s2){
  sum[0]=a[0]*s1+b[0]*s2;
  sum[1]=a[1]*s1+b[1]*s2;
  sum[2]=a[2]*s1+b[2]*s2;
  return sum;
}

/* structure for dihedral calculation */
typedef struct tagdihcalc_t{
  real phi; /* cis is zero, clockwise positive */
  real cosphi; /* cos(m, n) */
  real sign; /* (0, pi) is 1.0, otherwise -1.0 */

  real grad2;
  real grad[4][3]; /* gradient for each particle */

  real div; /* the divengent */
  real d4ij, d4ik, d4jj, d4jk, d4jl, d4kk, d4kl;

  unsigned int flags; /* a copy of flags used */
  int t1, t2, t3; /* gromacs shift indices */
  void *pbc; /* periodic boundary condition descriptor */
  int (*pbc_rv3_diff)(const void *, const real *xi, const real *xj, real *xij); /* a function to handle pbc */
    /* parameter order follows from the gromacs convention: the last is the difference
     * between the first two */
}dihcalc_t;

#define DIHCALC_GRAD  0x0001
#define DIHCALC_DIV   0x0002
/*#define DIHCALC_CONJ  0x0004 */
/*#define DIHCALC_PROJ  0x0008 */
#define DIHCALC_I     0x0010
#define DIHCALC_J     0x0020
#define DIHCALC_K     0x0040
#define DIHCALC_L     0x0080
#define DIHCALC_FOUR  (DIHCALC_I|DIHCALC_J|DIHCALC_K|DIHCALC_L)
/* the four atoms involved */
#define DIHCALC_ALL   (DIHCALC_FOUR|DIHCALC_GRAD|DIHCALC_DIV)
/* only I and L, so no divergence */
#define DIHCALC_ENDS  (DIHCALC_GRAD|DIHCALC_I|DIHCALC_L)

/* Calculates the dihedral angle, its gradient and the divegence
   of a conjugate field (to the gradient).

   For simplest calculation, dc can be NULL and flags=0,
   in this case only the dihedral angle is computed.

   The gradient and divergent are optional, depending on the flags,
   DIHCALC_GRAD and DIHCALC_DIV.
   To calculate the divergence, the force must be calculated.

   If you want to use a special routine for treating
   periodic boundary condition, assign a function pointer
   to pbc_rv3_diff before calling; and also supply additional
   information by the `pbc' entry.
   Otherwise, pbc_rv3_diff *must* be set to NULL.

   The procedures of calculating the angle and force
   are similar to that in GROMACS.

   We also introduce the conjugate field v
     v = grad(phi) / [grad(phi).grad(phi)],
   such that v.grad(phi) = 1.0.
   This field is not explicitly computed (since it's just a simple
   scaling of the gradient).
   The denominator is saved to grad2;

   You include only a few components of the gradient in calculating
   the conjugate field by passing `flags' a combination of
   DIHCALC_I, DIHCALC_J, DIHCALC_K, and DIHCALC_L.

   We also calculate the divergence of v, or div.v,
   It turns out the divergence of the gradient itself is always zero,
     div.grad(phi) = 0.
   Thus, we have
     div.v = -2 [ grad(phi).(grad grad(phi)).grad(phi) ] /[grad(phi).grad(phi)]^2.
   Components in grad(phi).(grad grad(phi)).grad(phi) involve the
   four gradient terms: d4ij, d4ik, ..., d4kl.

   Note, no matter what combination of DIHCALC_I, DIHCALC_J, DIHCALC_K, and DIHCALC_L
   is used, *all* moments d4ij, d4ik, ... , d4kl are calculated for safety.
   But only the involved ones are used to combined to produce the divergence.
 */
#define rv3_calcdihv(dc, x, idx, flags) rv3_calcdih(dc,x[idx[0]],x[idx[1]],x[idx[2]],x[idx[3]],flags)
ZCSTRCLS real rv3_calcdih(dihcalc_t *dc,
    const real *xi, const real *xj, const real *xk, const real *xl,
    unsigned int flags)
{
  real dot,scl,tol,vol,phi,sign,cosphi;
  real nxkj,nxkj2,m2,n2;
  real xij[3],xkj[3],xkl[3];
  real m[3], n[3]; /* the planar vector of xij x xkj,  and xkj x xkj */
  const real cosmax=(real)(1.0-6e-8);

  if(dc!=NULL && dc->pbc_rv3_diff != NULL){ /* handle pbc */
    dc->t1 = (*dc->pbc_rv3_diff)(dc->pbc, xi, xj, xij);
    dc->t2 = (*dc->pbc_rv3_diff)(dc->pbc, xk, xj, xkj);
    dc->t3 = (*dc->pbc_rv3_diff)(dc->pbc, xk, xl, xkl);
  }else{
    rv3_diff(xij, xi, xj);
    rv3_diff(xkj, xk, xj);
    rv3_diff(xkl, xk, xl);
  }
  nxkj2 = rv3_sqr(xkj);
  nxkj = (real)sqrt(nxkj2);
  tol=nxkj2*6e-8f;

  rv3_cross(m, xij, xkj);
  m2=rv3_sqr(m);
  rv3_cross(n, xkj, xkl);
  n2=rv3_sqr(n);
  if(m2>tol && n2>tol){
    scl=(real)sqrt(m2*n2);
    dot=rv3_dot(m, n);
    cosphi=dot/scl;
    if(cosphi>cosmax) cosphi=cosmax; else if(cosphi<-cosmax) cosphi=-cosmax;
  }else{
    cosphi=cosmax;
  }
  phi=(real)acos(cosphi);
  vol=rv3_dot(n, xij);
  sign=((vol>0.0f)?1.0f:(-1.0f));
  phi*=sign;
  if(dc != NULL){
    dc->phi=phi;
    dc->sign=sign;
    dc->cosphi=cosphi;
    dc->flags=flags;
  }

  /* optionally calculate the gradient */
  if(dc!=NULL && (flags&(DIHCALC_GRAD|DIHCALC_DIV)) ){ /* do gradient if divergence is required */
    /* clear divergence */
    dc->div=dc->d4ij=dc->d4ik=dc->d4jj=dc->d4jk=dc->d4jl=dc->d4kk=dc->d4kl=0.0f;

    /* calculate the gradient of the force
     * the direction of which increase the dihedral by one unit.
     */

    if(m2 > tol && n2 > tol){
      real gi[3],gj[3],gk[3],gl[3];
      real uvec[3],vvec[3],svec[3],p,q;
      real gi2, gj2, gk2, gl2, g2all,invg2;
      unsigned doi,doj,dok,dol;

      doi=(flags&DIHCALC_I);
      doj=(flags&DIHCALC_J);
      dok=(flags&DIHCALC_K);
      dol=(flags&DIHCALC_L);

      scl=nxkj/m2;
      rv3_smul2(gi, m, scl);
      scl=-nxkj/n2;
      rv3_smul2(gl, n, scl);

      p=rv3_dot(xij, xkj);
      p/=nxkj2;
      rv3_smul2(uvec, gi, p);
      q=rv3_dot(xkl, xkj);
      q/=nxkj2;
      rv3_smul2(vvec, gl, q);
      rv3_diff(svec, uvec, vvec);

      rv3_diff(gj, svec, gi);
      rv3_nsum2(gk, gl, svec);

      rv3_copy(dc->grad[0], gi);
      rv3_copy(dc->grad[1], gj);
      rv3_copy(dc->grad[2], gk);
      rv3_copy(dc->grad[3], gl);

      gi2=rv3_sqr(gi);
      gj2=rv3_sqr(gj);
      gk2=rv3_sqr(gk);
      gl2=rv3_sqr(gl);
      g2all=0.0f;
      if(doi) g2all+=gi2;
      if(doj) g2all+=gj2;
      if(dok) g2all+=gk2;
      if(dol) g2all+=gl2;
      dc->grad2=g2all;
      invg2=1.0f/g2all;

      if(flags&DIHCALC_DIV){
        real xkjv[3], nvv[3], mvv[3];
        real gjxij,gjmvv,gjxkl,gjnvv;
        real gkmvv,gknvv,gkxkl,gkxij;
        real kivkj,klvkj,ljvkj,ijvkj;
        real kikl, ijlj;
        real tmp1, tmp2;
        real sinmn;

        rv3_smul2(mvv, m, 1.0f/m2);
        rv3_smul2(nvv, n, 1.0f/n2);
        rv3_smul2(xkjv, xkj, 1.0f/nxkj);

        sinmn=vol*nxkj/(m2*n2);

        ijvkj=rv3_dot(xij, xkjv);
        kivkj=nxkj-ijvkj;
        klvkj=rv3_dot(xkl, xkjv);
        ljvkj=nxkj-klvkj;

        ijlj=ijvkj*ljvkj;
        kikl=kivkj*klvkj;

        gjxij=rv3_dot(gj, xij);
        gjxkl=rv3_dot(gj, xkl);
        gjmvv=rv3_dot(gj, mvv);
        gjnvv=rv3_dot(gj, nvv);
        gkxij=rv3_dot(gk, xij);
        gkxkl=rv3_dot(gk, xkl);
        gkmvv=rv3_dot(gk, mvv);
        gknvv=rv3_dot(gk, nvv);

        tmp1=nxkj2*sinmn;
        tmp2=tmp1/m2;
        dc->d4ij=kikl*tmp2;
        dc->d4ik=ijlj*tmp2;
        tmp2=tmp1/n2;
        dc->d4jl=kikl*tmp2;
        dc->d4kl=ijlj*tmp2;

        dc->d4jj=-(gjxij*gjmvv+gjxkl*gjnvv)/nxkj
                +2.0f*(kivkj*gjmvv-klvkj*gjnvv)*(-kikl*sinmn);

        dc->d4jk=(gjxij*gkmvv+gjxkl*gknvv)/nxkj
              +(-(gjmvv*ljvkj+gkmvv*klvkj)*(ijvkj*kivkj)
                +(gjnvv*ijvkj+gknvv*kivkj)*(ljvkj*klvkj) )*sinmn;

        dc->d4kk=-(gkxkl*gknvv+gkxij*gkmvv)/nxkj
                +2.0f*(ljvkj*gknvv-ijvkj*gkmvv)*(ijlj*sinmn);

        /* summarize */
        if((flags&DIHCALC_FOUR)==DIHCALC_FOUR){
          tmp1=dc->d4jj+dc->d4kk;
          tmp2=dc->d4ij+dc->d4ik+dc->d4jk+dc->d4jl+dc->d4kl;
        }else{
          tmp1=tmp2=0.0f;
          if(doj){ tmp1+=dc->d4jj; }
          if(dok){ tmp1+=dc->d4kk; }
          if(doi&&doj) tmp2+=dc->d4ij;
          if(doi&&dok) tmp2+=dc->d4ik;
          if(doj&&dok) tmp2+=dc->d4jk;
          if(doj&&dol) tmp2+=dc->d4jl;
          if(dok&&dol) tmp2+=dc->d4kl;
        }
        dc->div = -2.0f*(tmp1+2.0f*tmp2)*(invg2*invg2);
      } /* do divengence */

    }else{ /* clear the gradients */
      int j;
      for(j=0; j<4; j++) rv3_zero(dc->grad[j]);
    }
  }

  return phi;
}

#endif /* ZCOM_RV3__ */
#endif /* ZCOM_RV3 */


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

ZCINLINE void rv2_zero(real *x){x[0]=0; x[1]=0; }
ZCINLINE void rv2_copy(real *x, const real *src){x[0]=src[0]; x[1]=src[1]; }

ZCINLINE real rv2_sqr(const real *x){return x[0]*x[0]+x[1]*x[1];}
ZCINLINE real rv2_norm(const real *x){return (real)sqrt(x[0]*x[0]+x[1]*x[1]);}

ZCINLINE real *rv2_normalize(real *x){
  real r=rv2_norm(x);
  if(r>0.0){
    r=1.0/r;
    x[0]*=r;
    x[1]*=r;
  }
  return x;
}

ZCINLINE real rv2_dot(const real *x, const real *y){return x[0]*y[0]+x[1]*y[1];}

ZCINLINE real rv2_cross(const real *x, const real *y){
  return x[0]*y[1]-x[1]*y[0];
}

ZCINLINE real *rv2_neg(real *x){
  x[0]-=x[0];
  x[1]-=x[1];
  return x;
}

ZCINLINE real *rv2_neg2(real *nx, const real *x){
  nx[0]=-x[0];
  nx[1]=-x[1];
  return nx;
}

ZCINLINE real *rv2_inc(real *x, const real *dx){
  x[0]+=dx[0];
  x[1]+=dx[1];
  return x;
}
ZCINLINE real *rv2_dec(real *x, const real *dx){
  x[0]-=dx[0];
  x[1]-=dx[1];
  return x;
}
ZCINLINE real *rv2_sinc(real *x, real *dx, real s){
  x[0]+=s*dx[0];
  x[1]+=s*dx[1];
  return x;
}
ZCINLINE real *rv2_smul(real *x, real s){
  x[0]*=s;
  x[1]*=s;
  return x;
}
ZCINLINE real *rv2_smul2(real *y, const real *x, real s){
  y[0]=x[0]*s;
  y[1]=x[1]*s;
  return y;
}

/* for in-place difference use rv3_dec */
ZCINLINE real *rv2_diff(real *diff, const real *a, const real *b){
  diff[0]=a[0]-b[0];
  diff[1]=a[1]-b[1];
  return diff;
}

/* sum = a+b, for in-place addition use rv3_inc */
ZCINLINE real *rv2_sum2(real *sum, const real *a, const real *b){
  sum[0]=a[0]+b[0];
  sum[1]=a[1]+b[1];
  return sum;
}

/* sum = -a-b */
ZCINLINE real *rv2_nsum2(real *sum, const real *a, const real *b){
  sum[0]=-a[0]-b[0];
  sum[1]=-a[1]-b[1];
  return sum;
}

ZCINLINE real *rv2_lincomb2(real *sum, const real *a, const real *b, real s1, real s2){
  sum[0]=a[0]*s1+b[0]*s2;
  sum[1]=a[1]*s1+b[1]*s2;
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


#ifdef  ZCOM_ENDIAN
#ifndef ZCOM_ENDIAN__
#define ZCOM_ENDIAN__

#include <stdio.h>

/* Correct endianness from the current operating system to the desired one
 * The desired endianness is specified by tobig, while
 * the endianness of the current operating system is automatically computed.
 * A conversion takes place only if the two differ.
 * The enidan-corrected variable is saved in output, however,
 * for in-place conversion, pass NULL to output. */
ZCSTRCLS unsigned char *zcom_fix_endian(void *output, void *input, size_t len, int tobig)
{
  size_t i, ir;
  static int sysbig=-1;
  unsigned char *fixed=(unsigned char *)output;
  unsigned char *p=(unsigned char *)input;

  if(sysbig<0){ /* initial determine the machine's endianess */
    unsigned int feff=0xFEFF;
    unsigned char *s;
    s=(unsigned char *)(&feff);
    sysbig = ((*s==0xFF)?0:1);
    fprintf(stderr, "The current system is %s-endian.\n", (sysbig?"big":"little"));
  }

  if(tobig==sysbig){ /* no conversion is required */
    if(fixed==NULL){
      return p;
    }else{
      for(i=0; i<len; i++) fixed[i]=p[i];
      return fixed;
    }
  }else{

    if(fixed==NULL){ /* in-place conversion */
      for(i=0; i<(len/2); i++){
        char ch;
        ir=len-i-1;
        ch=p[i];
        p[i]=p[ir];
        p[ir]=ch;
      }
    }else{ /* out-of-place conversion */
      for(i=0; i<len; i++){
        for(i=0; i<len; i++)
          fixed[len-i-1]=p[i];
      }
    }
  }
  return p;
}


#endif /* ZCOM_ENDIAN__ */
#endif /* ZCOM_ENDIAN */



