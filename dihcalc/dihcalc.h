#ifndef DIHCALC_H__
#define DIHCALC_H__

/* structure for dihedral calculation */
typedef struct {
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
} dihcalc_t;

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

#define rv3_calcdihv(dc, x, idx, flags) \
  rv3_calcdih(dc, x[idx[0]], x[idx[1]], x[idx[2]], x[idx[3]], flags)
real rv3_calcdih(dihcalc_t *dc,
    const real *xi, const real *xj, const real *xk, const real *xl,
    unsigned int flags);

#endif

