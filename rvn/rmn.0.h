#ifndef RMN_H__
#define RMN_H__



#ifndef FMN_T
#define FMN_T fmn_t
typedef float fmn_t[D][DM];
#endif

#ifndef DMN_T
#define DMN_T dmn_t
typedef double dmn_t[D][DM];
#endif

#ifndef RMN_T
#define RMN_T rmn_t
/* only the second dimension needs to be DM, the first can be D */
typedef real rmn_t[D][DM];
#endif



#define rmn_print(r, nm, fmt, nl) rmn_fprint(stdout, r, nm, fmt, nl)

INLINE void rmn_fprint(FILE *fp, real r[][DM], const char *nm,
    const char *fmt, int nl)
{
  int i, j;
  if (nm) fprintf(fp, "%s:%c", nm, (nl ? '\n' : ' '));
  for (i = 0; i < D; i++) {
    for (j = 0; j < D; j++) {
      fprintf(fp, fmt, r[i][j], nl);
    }
    fprintf(fp, "%s", (nl ? "\n" : "; "));
  }
}



INLINE rvn_t *rmn_make(real x[][DM], ...)
{
  int i, j;
  va_list vl;

  va_start(vl, x);
  for (i = 0; i < D; i++)
    for (j = 0; j < D; j++)
      x[i][j] = (real) va_arg(vl, double);
  va_end(vl);
  return x;
}




INLINE rvn_t *rmn_zero(real x[][DM])
{
  int i;

  for (i = 0; i < D; i++) rvn_zero(x[i]);
  return x;
}



INLINE rvn_t *rmn_one(real x[][DM])
{
  int i;

  for (i = 0; i < D; i++) {
    rvn_zero(x[i]);
    x[i][i] = (real) 1;
  }
  return x;
}



/* generate a random orthonormal (unitary) DxD matrix */
INLINE rvn_t *rmn_randuni(real a[][DM])
{
  int i, j;
  real dot;

  rvn_randdir0(a[0]);
  for (i = 1; i < D; i++) {
    rvn_randdir0(a[i]);
    /* normalize against previous vectors */
    for (j = 0; j < i; j++) {
      dot = rvn_dot(a[i], a[j]);
      rvn_normalize( rvn_sinc(a[i], a[j], -dot) );
    }
  }
  return a;
}



#endif
