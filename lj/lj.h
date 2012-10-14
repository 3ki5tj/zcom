#include "def.h"
#include "util.h"
#include "hist.h"
#ifndef LJ_H__
#define LJ_H__

typedef struct {
  int i, j, in;
  real phi, psi, xi, dx[3], dr2;
} ljpair_t;


#define LJ_SWALLPAIRS 0x100 /* flag of usesw, save all (including out-of-range) pairs */

typedef struct {
  int d; /* dimension = 3 */
  int n; /* number of particles */
  int dof; /* degrees of freedom */
  real rho;
  real l, vol; /* side length and volume */
  real rc2, rc, rcdef; /* real / preferred rc */

  real * RESTRICT x; /* reduced unit (0, 1) */
  real * RESTRICT v, * RESTRICT f;
  real epot0, epot, epots; /* potential energy: pure, with tail correction and shifted potential energy */
  int iepot;  /* integer energy for square-well potential */
  real ekin, tkin, etot;
  real vir; /* virial */
  real epot_shift, epot_tail, p_tail;
  double t;

  int usesw; /* switched potential */
  real rs, a4, a5, a6, a7; /* parameters */
  ljpair_t *pr;
  int npr;
  real lap, f2, *gdg, *xdg;

  int usesq; /* square well potential */
  real ra, ra2, rb, rb2; /* -1 for (ra, rb) */

  hist_t *rdf; /* histogram for radial distribution function */
  int rdfnfr; /* number of frames in rdf */
} lj_t;

lj_t *lj_open(int n, int d, real rho, real rcdef);
void lj_close(lj_t *lj);
INLINE int lj_writepos(lj_t *lj, const real *x, const real *v, const char *fn);
#define LJ_LOADBOX 0x10
INLINE int lj_readpos(lj_t *lj, real *x, real *v, const char *fn, unsigned flags);
real lj_energy(lj_t *lj);
real lj_force(lj_t *lj);
INLINE void lj_vv(lj_t *lj, real dt);

#define lj_shiftcom(lj, v)    md_shiftcom(v, lj->n, lj->d)
#define lj_shiftang(lj, x, v) md_shiftang(x, v, lj->n, lj->d)

INLINE void lj_vscale(lj_t *lj, real tp, real ekt) 
 { md_vscale(lj->v, lj->n * lj->d, lj->dof, tp, ekt, &lj->ekin, &lj->tkin); }

INLINE void lj_vrescale(lj_t *lj, real tp, real thermdt) 
 { md_vrescale(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin); }

INLINE void lj_vrescalex(lj_t *lj, real tp, real thermdt)
 { md_vrescalex(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin); }

INLINE int lj_mcvrescale(lj_t *lj, real tp, real thermdt)
 { return md_mcvrescale(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin); }

/* compute pressure */
INLINE real lj_calcp(lj_t *lj, real tp)
  { return lj->rho * tp + lj->vir / (lj->d * lj->vol) + lj->p_tail; }

INLINE void lj_initsw(lj_t *lj, real rs);

/* initialize square well potential */
INLINE void lj_initsq(lj_t *lj, real ra, real rb)
{
  lj->ra = ra;
  lj->ra2 = ra * ra;
  lj->rb = rb;
  lj->rb2 = rb * rb;
  lj->usesq = 1;
  lj_energy(lj);
}

/* compute reference thermal dynamics variables using the equation of states
   return the average potential energy
   *P:  pressure
   *Ar: Helmholtz free energy (potential part)
   *Gr: Gibbs free energy (potential part)
   Reference:
   J. Karl Johnson et al. The Lennard-Jones equation of states revisited,
   Molecular Physics (1993) Vol. 78, No 3, 591-618 */
INLINE double lj_eos3d(double rho, double T, double *P, double *Ar, double *Gr)
{
  const double x1  =  0.8623085097507421;
  const double x2  =  2.976218765822098;
  const double x3  = -8.402230115796038;
  const double x4  =  0.1054136629203555;
  const double x5  = -0.8564583828174598;
  const double x6  =  1.582759470107601;
  const double x7  =  0.7639421948305453;
  const double x8  =  1.753173414312048;
  const double x9  =  2.798291772190376e+03;
  const double x10 = -4.8394220260857657e-02;
  const double x11 =  0.9963265197721935;
  const double x12 = -3.698000291272493e+01;
  const double x13 =  2.084012299434647e+01;
  const double x14 =  8.305402124717285e+01;
  const double x15 = -9.574799715203068e+02;
  const double x16 = -1.477746229234994e+02;
  const double x17 =  6.398607852471505e+01;
  const double x18 =  1.603993673294834e+01;
  const double x19 =  6.805916615864377e+01;
  const double x20 = -2.791293578795945e+03;
  const double x21 = -6.245128304568454;
  const double x22 = -8.116836104958410e+03;
  const double x23 =  1.488735559561229e+01;
  const double x24 = -1.059346754655084e+04;
  const double x25 = -1.131607632802822e+02;
  const double x26 = -8.867771540418822e+03;
  const double x27 = -3.986982844450543e+01;
  const double x28 = -4.689270299917261e+03;
  const double x29 =  2.593535277438717e+02;
  const double x30 = -2.694523589434903e+03;
  const double x31 = -7.218487631550215e+02;
  const double x32 =  1.721802063863269e+02;
  const double gamma = 3.0;
  double a[8], b[6], c[8], d[6], G[6], F, rhop, rho2 = rho*rho, Pa = 0., Pb = 0., U;
  int i;
    
  a[0] = x1*T + x2*sqrt(T) + x3 + x4/T + x5/(T*T);
  a[1] = x6*T + x7 + x8/T + x9/(T*T);
  a[2] = x10*T + x11 + x12/T;
  a[3] = x13;
  a[4] = x14/T + x15/(T*T);
  a[5] = x16/T;
  a[6] = x17/T + x18/(T*T);
  a[7] = x19/(T*T);
  b[0] = (x20 + x21/T)/(T*T);
  b[1] = (x22 + x23/(T*T))/(T*T);
  b[2] = (x24 + x25/T)/(T*T);
  b[3] = (x26 + x27/(T*T))/(T*T);
  b[4] = (x28 + x29/T)/(T*T);
  b[5] = (x30 + x31/T + x32/(T*T))/(T*T);
  c[0] = x2*sqrt(T)/2 + x3 + 2*x4/T + 3*x5/(T*T);
  c[1] = x7 + 2*x8/T + 3*x9/(T*T);
  c[2] = x11 + 2*x12/T;
  c[3] = x13;
  c[4] = 2*x14/T + 3*x15/(T*T);
  c[5] = 2*x16/T;
  c[6] = 2*x17/T + 3*x18/(T*T);
  c[7] = 3*x19/(T*T);
  d[0] = (3*x20 + 4*x21/T)/(T*T);
  d[1] = (3*x22 + 5*x23/(T*T))/(T*T);
  d[2] = (3*x24 + 4*x25/T)/(T*T);
  d[3] = (3*x26 + 5*x27/(T*T))/(T*T);
  d[4] = (3*x28 + 4*x29/T)/(T*T);
  d[5] = (3*x30 + 4*x31/T + 5*x32/(T*T))/(T*T);
  
  F = exp(-gamma*rho*rho);
  G[0] = (1 - F)/(2*gamma);
  for (rhop = 1, i = 1; i < 6; i++) {
    rhop *= rho*rho;
    G[i] = -(F*rhop - 2*i*G[i-1])/(2*gamma);
  }

  if (Ar) *Ar = 0.;
  if (P)  Pa = Pb = 0;
  for (U = 0, i = 7; i >= 0; i--) {
    U = rho * (c[i]/(i+1) + U);
    if (Ar) *Ar = rho * (a[i]/(i+1) + (*Ar));
    if (P)  Pa  = rho * (a[i] + Pa);
  }
  
  for (i = 5; i >= 0; i--) {
    U += d[i]*G[i];
    if (Ar) *Ar += b[i]*G[i];
    if (P) Pb = rho2*(b[i] + Pb);
  }
  if (P) *P = rho*(T + Pa + F*Pb);
  if (Gr) *Gr = *Ar + *P/rho - T;
  return U;
}

#endif

