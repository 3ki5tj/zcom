#include "def.h"
#include "util.h"
#include "hist.h"
#include "md.h"
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

  int usesw; /* switched potential */
  real rs, a4, a5, a6, a7; /* parameters */
  ljpair_t *pr;
  int npr;
  real lap, f2, *gdg, *xdg;

  int usesq; /* square well potential */
  int esqinf;
  real ra, ra2, rb, rb2; /* -1 for (ra, rb) */
  real rmin; /* minimal pair distance */

  unsigned isclone; /* is a clone copy, don't free pointers */
} lj_t;

INLINE real lj_energy(lj_t *lj);
INLINE real lj_force(lj_t *lj);

#define lj_shiftcom(lj, v)    md_shiftcom(v, lj->n, lj->d)
#define lj_shiftang(lj, x, v) md_shiftang(x, v, lj->n, lj->d)



/* initialize a fcc lattice */
static void lj_initfcc2d(lj_t *lj)
{
  int i, j, id, n1, n = lj->n;
  real a;

  n1 = (int) (pow(2*n, 1.0/lj->d) + .999999); /* # of particles per side */
  a = 1.f/n1;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++) {
      if ((i+j) % 2 != 0) continue;
      lj->x[id*2 + 0] = (i + .5f) * a;
      lj->x[id*2 + 1] = (j + .5f) * a;
      id++;
    }
}



/* initialize a fcc lattice */
static void lj_initfcc3d(lj_t *lj)
{
  int i, j, k, id, n1, n = lj->n;
  real a;

  n1 = (int) (pow(2*n, 1.0/lj->d) + .999999); /* # of particles per side */
  a = 1.f/n1;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++)
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 != 0) continue;
        lj->x[id*3 + 0] = (i + .5f) * a;
        lj->x[id*3 + 1] = (j + .5f) * a;
        lj->x[id*3 + 2] = (k + .5f) * a;
        id++;
      }
}



/* set density and compute tail corrections */
INLINE void lj_setrho(lj_t *lj, real rho)
{
  double irc, irc3, irc6;
  int i;

  lj->rho = rho;
  lj->l = (real) pow(1.*lj->n/rho, 1./lj->d);
  for (lj->vol = 1.f, i = 0; i < lj->d; i++) lj->vol *= lj->l;
  if ((lj->rc = lj->rcdef) > lj->l * .5f) lj->rc = lj->l * .5f;
  lj->rc2 = lj->rc * lj->rc;
  irc = 1.0/lj->rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  if (lj->usesw) { /* assume u(L/2) = 0 */
    lj->epot_shift = 0.f;
    lj->epot_tail = 0.f;
    lj->p_tail = 0.f;
  } else {
    lj->epot_shift = (real)( 4*irc6*(irc6 - 1) );
    if (lj->d == 3) {
      lj->epot_tail = (real)( 8*M_PI*rho*lj->n/9*(irc6 - 3)*irc3 );
      lj->p_tail = (real)( 32*M_PI*rho*rho/9*(irc6 - 1.5)*irc3 );
    } else if (lj->d == 2) {
      lj->epot_tail = (real) (M_PI*rho*lj->n*(.4*irc6 - 1)*irc3*irc);
      lj->p_tail = (real) (M_PI*rho*rho*(1.6*irc6 - 2)*irc3*irc);
    }
  }
}



INLINE real lj_pbc(real x, real l)
  { return (real)((x - ((int)((x)+1000.5) - 1000))*l) ; }



INLINE real *lj_vpbc2d(real *v, real l)
  { v[0] = lj_pbc(v[0], l); v[1] = lj_pbc(v[1], l); return v; }



INLINE real *lj_vpbc3d(real *v, real l)
  { v[0] = lj_pbc(v[0], l); v[1] = lj_pbc(v[1], l); v[2] = lj_pbc(v[2], l); return v; }



INLINE real lj_pbcdist2_3d(real *dx, const real *a, const real *b, real l)
  { return rv3_sqr(lj_vpbc3d(rv3_diff(dx, a, b), l)); }



/* create an open structure */
lj_t *lj_open(int n, int d, real rho, real rcdef)
{
  lj_t *lj;
  int i;

  xnew(lj, 1);
  lj->n = n;
  lj->d = d;
  lj->dof = n * d - d * (d+1)/2;
  xnew(lj->f, n * d);
  xnew(lj->v, n * d);
  xnew(lj->x, n * d);

  lj->rcdef = rcdef;
  lj_setrho(lj, rho);

  lj->esqinf = 1000000;

  if (lj->d == 3) lj_initfcc3d(lj); else lj_initfcc2d(lj);

  /* init. random velocities */
  for (i = 0; i < n * d; i++) lj->v[i] = (real) (rnd0() - .5);

  lj_shiftcom(lj, lj->v);
  lj_shiftang(lj, lj->x, lj->v);
  lj->ekin = md_ekin(lj->v, lj->n * lj->d, lj->dof, &lj->tkin);

  lj->isclone = 0;

  lj_force(lj);
  return lj;
}



/* copy flags */
#define LJ_CPX   0x0001
#define LJ_CPV   0x0002
#define LJ_CPF   0x0004
#define LJ_CPPR  0x0020
#define LJ_CPGDG 0x0040
#define LJ_CPXDG 0x0080
#define LJ_CPXVF (LJ_CPX|LJ_CPV|LJ_CPF)

#define lj_copyvec(lj, t, s) memcpy(t, s, lj->d * lj->n * sizeof(real))

/* copy from src to dest
 * cannot copy vectors other than xvf */
INLINE lj_t *lj_copy(lj_t *dest, const lj_t *src, unsigned flags)
{
  /* to preserve the pointers before the memcpy(dest, src) call */
  real *x = dest->x, *v = dest->v, *f = dest->f;

  memcpy(dest, src, sizeof(lj_t));

  if (flags & LJ_CPX) lj_copyvec(src, x, src->x);
  dest->x = x;
  if (flags & LJ_CPV) lj_copyvec(src, v, src->v);
  dest->v = v;
  if (flags & LJ_CPF) lj_copyvec(src, f, src->f);
  dest->f = f;
  return dest;
}



/* make new copy */
INLINE lj_t *lj_clone(const lj_t *src, unsigned flags)
{
  int nd = src->n * src->d;
  lj_t *dest;

  xnew(dest, 1);
  memcpy(dest, src, sizeof(lj_t));
  /* unless specified in flags,
   * arrays are copied literally as pointers */
  dest->isclone = LJ_CPPR | LJ_CPGDG | LJ_CPXDG;
  if (flags & LJ_CPX) {
    xnew(dest->x, nd);
    lj_copyvec(src, dest->x, src->x);
  } else {
    dest->isclone |= LJ_CPX;
  }
  if (flags & LJ_CPV) {
    xnew(dest->v, nd);
    lj_copyvec(src, dest->v, src->v);
  } else {
    dest->isclone |= LJ_CPV;
  }
  if (flags & LJ_CPF) {
    xnew(dest->f, nd);
    lj_copyvec(src, dest->f, src->v);
  } else {
    dest->isclone |= LJ_CPF;
  }
  return dest;
}



void lj_close(lj_t *lj)
{
  if ( !(lj->isclone & LJ_CPX) ) free(lj->x);
  if ( !(lj->isclone & LJ_CPV) ) free(lj->v);
  if ( !(lj->isclone & LJ_CPF) ) free(lj->f);
  if ( !(lj->isclone & LJ_CPPR)  && lj->pr)
    free(lj->pr);
  if ( !(lj->isclone & LJ_CPGDG) && lj->gdg)
    free(lj->gdg);
  if ( !(lj->isclone & LJ_CPXDG) && lj->xdg)
    free(lj->xdg);
  free(lj);
}



/* write position (and velocity)
 * Note 1: *actual* position, not unit position is used
 * Note 2: coordinates are *not* wrapped back into the box */
INLINE int lj_writepos(lj_t *lj, const real *x, const real *v, const char *fn)
{
  FILE *fp;

  if (fn == NULL) fn = "lj.pos";
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d %d %.14e\n", lj->d, lj->n, (v != NULL), lj->l);
  md_writepos(fp, x, v, lj->n, lj->d, lj->l);
  fclose(fp);
  return 0;
}



#define LJ_LOADBOX 0x10

/* read position file (which may include velocity) */
INLINE int lj_readpos(lj_t *lj, real *x, real *v, const char *fn, unsigned flags)
{
  char s[1024];
  FILE *fp;
  int i, j, ret = -1, hasv = 0;
  double l0;

  if (fn == NULL) fn = "lj.pos";
  xfopen(fp, fn, "r", return -1);

  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') { /* simplified format, old version */
    fprintf(stderr, "Warning: %s has no information line\n", fn);
    rewind(fp);
  } else {
    if (4 != sscanf(s + 1, "%d%d%d%lf", &i, &j, &hasv, &l0)
        || i != lj->d || j != lj->n) {
      fprintf(stderr, "first line is corrupted:\n%s", s);
      goto ERR;
    }
    if (fabs(l0 - lj->l) > 1e-5*lj->l) { /* verify the box size */
      if (flags & LJ_LOADBOX) {
        lj->l = (real) l0;
        for (lj->vol = 1, j = 0; j < lj->d; j++) lj->vol *= lj->l;
        lj_setrho(lj, lj->n/lj->vol);
      } else {
        fprintf(stderr, "box mismatch l %g, should be %g\n", l0, lj->l);
        goto ERR;
      }
    }
  }

  ret = md_readpos(fp, x, hasv ? v : NULL, lj->n, lj->d, lj->l);
ERR:
  fclose(fp);
  return ret;
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



/* initialize square well potential */
INLINE void lj_initsq(lj_t *lj, real ra, real rb)
{
  lj->ra2 = (lj->ra = ra) * ra;
  lj->rb2 = (lj->rb = rb) * rb;
  lj->usesq = 1;
  lj_energy(lj);
}



/* initialize coefficients for the switched potential */
INLINE void lj_initsw(lj_t *lj, real rs)
{
  real rs2, rs3, rs6, rs15, dr, dr2, dr3, dr4, f1, f2, f26, f13;

  lj->rs = rs;
  dr = lj->rs - lj->rc;
  die_if (dr > 0, "rs %g, rc %g\n", lj->rs, lj->rc);

  rs2 = rs*rs;
  rs3 = rs2*rs;
  rs6 = rs3*rs3;
  rs15 = rs6*rs6*rs3;
  dr2 = dr*dr;
  dr3 = dr2*dr;
  dr4 = dr3*dr;
  f1 = rs6 - 1.f;
  f2 = rs6 - 2.f;
  f13 = 2.f*rs6 - 13.f;
  f26 = 7.f*rs6 - 26.f;

  f1 *= rs3;
  f2 *= dr*rs2;
  f13 *= dr3;
  f26 *= dr2*rs;
  lj->a4 = -4.f*(35.f*f1 + 90.f*f2 + 28.f*f13 + 15.f*f26)/(dr4*rs15);
  lj->a5 = 24.f*(14.f*f1 + 39.f*f2 + 14.f*f13 + 7.f*f26)/(dr2*dr3*rs15);
  lj->a6 = -4.f*(70.f*f1 + 204.f*f2 + 84.f*f13 + 39.f*f26)/(dr3*dr3*rs15);
  lj->a7 = 16.f*(5.f*f1 + 15.f*f2 + 7.f*f13 + 3.f*f26)/(dr3*dr4*rs15);

  xnew(lj->pr, lj->n * lj->n);
  xnew(lj->gdg, lj->n * lj->d);
  lj->npr = 0;
  lj->usesw = 1;
}



/* compute the switch potential phi(r) and its derivatives
 * fscal = -phi = uij'/rij
 * psi = phi'/rij
 * xi = psi'/rij
 * Laplacian = psi*rij^2 + 3*phi = psi*rij^2 - 3*fscal */
INLINE real lj_potsw(lj_t *lj, real r, real *fscal, real *psi, real *xi)
{
  if (r < lj->rs) { /* normal lj */
    real invr2, invr6, invr8;
    invr2 = 1 / (r*r);
    invr6 = invr2 * invr2 * invr2;
    invr8 = invr6 * invr2;
    *fscal = (48 * invr6 - 24) * invr8;
    *psi = (672 * invr6 - 192) * invr8 * invr2;
    *xi = -(10752 * invr6 - 1920) * invr6 * invr6;
    return 4 * invr6 * (invr6 - 1);
  } else if (r < lj->rc) { /* polynomial */
    real dr, dr2, dr3, fs, ps, xs, invr, invr2;
    real a4 = lj->a4, a5 = lj->a5, a6 = lj->a6, a7 = lj->a7;
    invr = 1/r;
    dr = r - lj->rc;
    invr2 = invr * invr;
    dr2 = dr * dr;
    dr3 = dr2 * dr;
    fs = dr3*(4*a4 + dr*(5*a5 + dr*(6*a6 + dr*7*a7)))*invr;
    *fscal = -fs;
    ps = dr2*(12*a4 + dr*(20*a5 + dr*(30*a6 + dr*42*a7)));
    *psi = (ps - fs)*invr2;
    xs = dr*(24*a4 + dr*(60*a5 + dr*(120*a6 + dr*210*a7)));
    *xi = (xs*invr - 3*(*psi))*invr2;
    return (dr2*dr2)*(a4 + dr*(a5 + dr*(a6 + dr*a7)));
  } else { /* out of range */
    *fscal = 0;
    *psi = 0;
    *xi = 0;
    return 0;
  }
}



/* 3D energy for square, lj members are not altered */
static int lj_energysq3d(lj_t *lj, rv3_t *x, real *rmin)
{
  real dx[3], dr2, ra2 = lj->ra2, rb2 = lj->rb2, l = lj->l, rm2 = 1e30;
  int i, j, iu = 0, n = lj->n, col = 0;

  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], l);
      if (dr2 < ra2) {
        iu += lj->esqinf;
        col++;
      } else if (dr2 < rb2) {
        iu--;
      }
      if (dr2 < rm2) rm2 = dr2;
    }
  }
  if (rmin != NULL) *rmin = (real) sqrt(rm2);
  /* force the energy to zero in the hard sphere case */
  if (fabs(ra2 - rb2) < 1e-6 && col == 0) {
    iu = 0;
  }
  return iu;
}



/* compute 3D energy for switched potential */
INLINE real lj_energysw3d(lj_t *lj, rv3_t *x, real *virial, real *laplace)
{
  int i, j, n = lj->n;
  real dx[3], dr2, dr, l = lj->l, d = (real) lj->d;
  real fscal, psi, xi, ep, vir, lap;

  ep = lap = vir = 0.f;
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], l);
      if (dr2 > lj->rc2) continue;

      dr = (real) sqrt(dr2);
      ep += lj_potsw(lj, dr, &fscal, &psi, &xi);
      lap += 2.f*(psi*dr2 - d*fscal);
      vir += fscal*dr2;
      rv3_smul(dx, fscal);
    }
  }
  if (virial) *virial = vir;
  if (laplace) *laplace = lap;
  return ep;
}



/* 3D compute force and virial, return energy */
static real lj_energylj3d(lj_t *lj, rv3_t *x, real *virial, real *ep0, real *eps)
{
  real dx[3], dr2, dr6, ep, vir, l = lj->l, rc2 = lj->rc2;
  int i, j, prcnt = 0, n = lj->n;

  if (virial) *virial = 0.f;
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], l);
      if (dr2 > rc2) continue;
      dr2 = 1.f / dr2;
      dr6 = dr2 * dr2 * dr2;
      vir += dr6 * (48.f * dr6 - 24.f); /* f.r */
      ep += 4 * dr6 * (dr6 - 1);
      prcnt++;
    }
  }
  if (ep0) *ep0 = ep;
  if (eps) *eps = ep - prcnt * lj->epot_shift; /* shifted energy */
  if (virial) *virial = vir;
  return ep + lj->epot_tail; /* unshifted energy */
}



/* energy evaluation, do not change members of `lj' */
INLINE real lj_energyx3d(lj_t *lj, rv3_t *x, real *vir, int *iep, real *rmin,
    real *ep0, real *eps, real *lap)
{
  real u;
  if (lj->usesq) {
    *iep = lj_energysq3d(lj, x, rmin);
    u = (real) (*iep);
  } else if (lj->usesw) {
    u = lj_energysw3d(lj, x, vir, lap);
  } else {
    u = lj_energylj3d(lj, x, vir, ep0, eps);
  }
  if (lj->usesq || lj->usesw) {
    if (ep0) *ep0 = u;
    if (eps) *eps = u;
  }
  return u;
}



/* energy evaluation, do not change members of `lj' */
INLINE real lj_energyx(lj_t *lj, real *x, real *vir, int *iep, real *rmin,
    real *ep0, real *eps, real *lap)
{
  return lj->d == 2 ?
      lj_energyx2d(lj, (rv2_t *) x, vir, iep, rmin, ep0, eps, lap) :
      lj_energyx3d(lj, (rv3_t *) x, vir, iep, rmin, ep0, eps, lap);
}



/* compute the energy of the current configuration and set lj->epot */
INLINE real lj_energy3d(lj_t *lj)
{
  return lj->epot = lj_energyx3d(lj, (rv3_t *) lj->x, &lj->vir, &lj->iepot,
      &lj->rmin, &lj->epot0, &lj->epots, &lj->lap);
}



INLINE real lj_energy(lj_t *lj)
{
  return (lj->d == 2) ? lj_energy2d(lj) : lj_energy3d(lj);
}



/* compute 3D switched force, save derivative information to lj->pr */
INLINE real lj_forcesw3d(lj_t *lj, rv3_t *x, rv3_t *f, ljpair_t *pr,
    int *ljnpr, real *virial, real *f2, real *laplace)
{
  int i, j, n = lj->n, npr;
  real dx[3], dr2, dr, l = lj->l, d = (real) lj->d;
  real fscal, psi, xi, ep, vir, lap;

  npr = 0;
  ep = lap = vir = 0.f;
  for (i = 0; i < n; i++) rv3_zero(f[i]);
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], l);
      if (dr2 > lj->rc2) {
        /* save out-of-range pairs, so we can later
         * locate the pair from indices i and j using getpairindex() */
        if (lj->usesw & LJ_SWALLPAIRS) {
          rv3_copy(pr->dx, dx);
          pr->i = i;
          pr->j = j;
          pr->phi = pr->psi = pr->xi = 0.f;
          pr->dr2 = dr2;
          pr->in = 0;
          pr++; npr++;
        }
        continue;
      }

      rv3_copy(pr->dx, dx);
      pr->dr2 = dr2;
      dr = (real) sqrt(dr2);
      ep += lj_potsw(lj, dr, &fscal, &psi, &xi);
      pr->phi = -fscal; /* phi = u'/r */
      pr->psi = psi;  /* psi = phi'/r */
      pr->xi = xi;  /* xi = psi'/r */
      lap += 2.f*(psi*dr2 - d*fscal);
      vir += fscal*dr2; /* f.r */
      rv3_smul(dx, fscal);
      pr->i = i;
      pr->j = j;
      rv3_inc(f[i], dx);
      rv3_dec(f[j], dx);
      pr->in = 1;
      pr++; npr++;
    }
  }
  if (ljnpr) *ljnpr = npr;
  if (virial) *virial = vir;
  if (laplace) *laplace = lap;
  if (f2) for (*f2 = 0.0, i = 0; i < n; i++) *f2 += rv3_sqr(f[i]);
  return ep;
}



/* 3D compute force and virial, return energy */
static real lj_forcelj3d(lj_t *lj, rv3_t *x, rv3_t *f, real *virial,
    real *ep0, real *eps, real *f2, real *laplace)
{
  real dx[3], fi[3], dr2, dr6, fs, tmp, ep, vir, lap, l = lj->l, rc2 = lj->rc2;
  int i, j, d, prcnt = 0, n = lj->n;

  for (i = 0; i < n; i++) rv3_zero(f[i]);
  for (ep = vir = lap = 0, i = 0; i < n - 1; i++) {
    rv3_zero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], l);
      if (dr2 > rc2) continue;
      dr2 = 1.f / dr2;
      dr6 = dr2 * dr2 * dr2;
      fs = dr6 * (48.f * dr6 - 24.f); /* f.r */
      vir += fs; /* f.r */
      if (laplace) /* 2.f for it applies to both particles */
        lap += 2.f * ((168 - 12*3) * dr6 - (48 - 6*3)) * dr6 * dr2;

      fs *= dr2; /* f.r / r^2 */
      for (d = 0; d < 3; d++) {
        tmp = dx[d] * fs;
        fi[d] += tmp;
        f[j][d] -= tmp;
      }

      ep += 4 * dr6 * (dr6 - 1);
      prcnt++;
    }
    rv3_inc(f[i], fi);
  }
  if (ep0) *ep0 = ep;
  if (eps) *eps = ep - prcnt * lj->epot_shift; /* shifted energy */
  if (virial) *virial = vir;
  if (laplace) *laplace = 4*lap;
  if (f2) for (*f2 = 0.f, i = 0; i < n; i++) *f2 += rv3_sqr(f[i]);
  return ep + lj->epot_tail; /* unshifted energy */
}



INLINE real lj_force3d(lj_t *lj)
{
  if (lj->usesq) return lj->epot = lj->epot0 = lj->epots = (real) lj_energysq3d(lj,
    (rv3_t *) lj->x, &lj->rmin); /* no force for square well */
  if (lj->usesw) return lj->epot = lj->epot0 = lj->epots = lj_forcesw3d(lj,
    (rv3_t *) lj->x, (rv3_t *) lj->f,
    lj->pr, &lj->npr, &lj->vir, &lj->f2, &lj->lap);
  return lj->epot = lj_forcelj3d(lj, (rv3_t *) lj->x, (rv3_t *) lj->f,
    &lj->vir, &lj->epot0, &lj->epots, &lj->f2, &lj->lap);
}



INLINE real lj_force(lj_t *lj)
{
  return (lj->d == 2) ? lj_force2d(lj) : lj_force3d(lj);
}



/* calculate the configurational temperature (bc) for switched potential
 * bc = div(v), where v = g/(g.g), g = grad U,
 * udb = v . grad bc
 * bvir = x . grad bc */
INLINE real lj_bconfsw3d(lj_t *lj, real *udb)
{
  int i, j, ipr, npr = lj->npr, n = lj->n;
  ljpair_t *pr;
  real dg[3], dh[3];
  real phi, psi, xi, d = (real) lj->d;
  real dgdx, dg2, dr2, m = 0.f, h2;
  real gdlap = 0.f, gdm = 0.f, bc, invg2, invg4;
  real dlap, dgdx2;
  rv3_t *h = (rv3_t *) lj->gdg, *f = (rv3_t *) lj->f;

  if (udb) for (i = 0; i < n; i++) rv3_zero(h[i]);

  for (ipr = 0; ipr < npr; ipr++) {
    pr = lj->pr + ipr;
    i = pr->i;
    j = pr->j;
    phi = pr->phi;
    psi = pr->psi;

    dg2 = rv3_sqr( rv3_diff(dg, f[j], f[i]) );
    dgdx = rv3_dot(dg, pr->dx);
    m += psi*(dgdx2 = dgdx*dgdx) + phi*dg2; /* M = g g : grad grad U */
    if (udb) {
      dr2 = pr->dr2;
      xi = pr->xi;
      dlap = xi*dr2 + (2.f + d)*psi;
      gdlap += dlap * dgdx; /* 0.5 g . grad laplace U */
      gdm += (xi*dgdx2 + 3.f*psi*dg2)*dgdx; /* g g g : grad grad grad U, first larger */
      rv3_lincomb2(dh, pr->dx, dg, psi * dgdx, phi);
      rv3_sinc(h[i], dh,  2.f);
      rv3_sinc(h[j], dh, -2.f);
    }
  }
  m *= 2.f;
  gdlap *= 2.f;
  gdm *= 2.f;
  invg2 = 1.f/lj->f2;
  invg4 = invg2 * invg2;
  bc = (lj->lap - m*invg2)*invg2; /* configuration temperature */

  if (udb) {
    for (h2 = 0.f, i = 0; i < n; i++) h2 += rv3_sqr(h[i]);
    /* (1/g) \partial^2 g/ \partial E^2 = <bc*bc + udb>
       \partial bc/ \partial E = <d(bc)^2 + udb> */
    *udb = invg4*(gdlap - (lj->lap*m + h2 + gdm)*invg2 + 2.f*m*m*invg4);
  }

  return bc;
}



/* return r r : grad grad U, must be called after force */
INLINE real lj_vir2sw3d(lj_t *lj)
{
  int ipr, npr = lj->npr;
  real vir2 = 0.f;

  for (ipr = 0; ipr < npr; ipr++) {
    ljpair_t *pr = lj->pr + ipr;
    vir2 += (pr->psi * pr->dr2 + pr->phi) * pr->dr2;
  }
  return vir2;
}



/* compute pressure */
INLINE real lj_calcp(lj_t *lj, real tp)
{ return (lj->dof * tp + lj->vir) / (lj->d * lj->vol) + lj->p_tail; }



/* compute pressure, ideal gas part from the kinetic energy  */
INLINE real lj_calcpk(lj_t *lj)
{ return (2.f * lj->ekin + lj->vir) / (lj->d * lj->vol) + lj->p_tail; }



#endif

