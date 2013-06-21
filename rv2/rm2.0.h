#ifndef RM2_H__
#define RM2_H__



#ifndef FM2_T
#define FM2_T fm2_t
typedef float fm2_t[2];
#endif

#ifndef DM2_T
#define DM2_T dm2_t
typedef double dm2_t[2];
#endif

#ifndef RM2_T
#define RM2_T rm2_t
typedef real rm2_t[2];
#endif



#define rm2_print(r, nm, fmt, nl) rm2_fprint(stdout, r, nm, fmt, nl)

INLINE void rm2_fprint(FILE *fp, real r[2][2], const char *nm,
    const char *fmt, int nl)
{
  int i, j;
  if (nm) fprintf(fp, "%s:%c", nm, (nl ? '\n' : ' '));
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      fprintf(fp, fmt, r[i][j], nl);
    }
    fprintf(fp, "%s", (nl ? "\n" : "; "));
  }
}



INLINE rv2_t *rm2_make(real x[2][2], real a00, real a01, real a10, real a11)
{
  rv2_make(x[0], a00, a01);
  rv2_make(x[1], a10, a11);
  return x;
}



#define rm2_makem(rx, x) rm2_make(rx, (real) x[0][0], (real) x[0][1], \
    (real) x[1][0], (real) x[1][1])



#define rm2_zero(x) rm2_make(x, 0, 0, 0, 0)



#define rm2_one(x) rm2_make(x, 1, 0, 0, 1)



/* generate a random orthonormal (unitary) 2x2 matrix */
INLINE rv2_t *rm2_rnduni(real a[2][2])
{
  rv2_rnddir0(a[0]);
  rv2_make(a[1], a[0][1], -a[0][0]);
  if (rnd0() > 0.5) rv2_neg(a[1]);
  return a;
}



#endif
