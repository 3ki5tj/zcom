#include "md.c"
#ifndef LJ_C__
#define LJ_C__
#include "lj.h"

/* set density and compute tail corrections */
static void lj_setrho(lj_t *lj, real rho)
{
  double irc, irc3, irc6;
  int i;

  lj->rho = rho;
  lj->l = (real) pow(1.*lj->n/rho, 1./lj->d);
  for (lj->vol = 1.f, i = 0; i < lj->d; i++) lj->vol *= lj->l;
  if ((lj->rc = lj->rcdef) > lj->l * .5f) lj->rc = lj->l * .5f;
  lj->rc2 = lj->rc * lj->rc;
  irc = 1.0/lj->rc;
  irc3 = irc*irc*irc; irc6 = irc3*irc3;
  if (lj->usesw) { /* assume u(L/2) = 0 */
    lj->epot_shift = 0.f;
    lj->epot_tail = 0.f;
    lj->p_tail = 0.f;
  } else {
    lj->epot_shift = (real)( 4*irc6*(irc6-1) );
    if (lj->d == 3) {
      lj->epot_tail = (real)( 8*M_PI*rho*lj->n/9*(irc6 - 3)*irc3 );
      lj->p_tail = (real)( 32*M_PI*rho*rho/9*(irc6 - 1.5)*irc3 );
    } else if (lj->d == 2) {
      lj->epot_tail = (real) (M_PI*rho*lj->n*(.4*irc6 - 1)*irc3*irc);
      lj->p_tail = (real) (M_PI*rho*rho*(1.6*irc6 - 2)*irc3*irc);
    }
  }
}

/* initialize coefficients for the switched potential */
void lj_initsw(lj_t *lj, real rs)
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
 * laplacian = psi*rij^2 + 3*phi = psi*rij^2 - 3*fscal */
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

INLINE real lj_pbc(real x, real l)
  { return (real)((x - ((int)((x)+1000.5) - 1000))*l) ; }

INLINE real *lj_vpbc2d(real *v, real l)
  { v[0] = lj_pbc(v[0], l); v[1] = lj_pbc(v[1], l); return v; }

INLINE real *lj_vpbc3d(real *v, real l)
  { v[0] = lj_pbc(v[0], l); v[1] = lj_pbc(v[1], l); v[2] = lj_pbc(v[2], l); return v; }

INLINE real lj_pbcdist2_3d(real *dx, const real *a, const real *b, real l)
  { return rv3_sqr(lj_vpbc3d(rv3_diff(dx, a, b), l)); }

/* 3D energy for square */
static int lj_energysq3d(lj_t *lj, rv3_t *x)
{
  real dx[3], dr2, ra2 = lj->ra2, rb2 = lj->rb2, l = lj->l;
  int i, j, iu = 0, n = lj->n;

  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], l);
      if (dr2 < ra2)
        iu += 1000000;
      else if (dr2 < rb2)
        iu--;
    }
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
  real dx[3], dr2, dr6, ep, l = lj->l, rc2 = lj->rc2;
  int i, j, prcnt = 0, n = lj->n;

  if (virial) *virial = 0.f;
  for (ep = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], l);
      if (dr2 > rc2) continue;
      dr2 = 1.f / dr2;
      dr6 = dr2 * dr2 * dr2;
      if (virial) *virial += dr6 * (48.f * dr6 - 24.f); /* f.r */
      ep += 4 * dr6 * (dr6 - 1);
      prcnt++;
    }
  }
  if (ep0) *ep0 = ep;
  if (eps) *eps = ep - prcnt * lj->epot_shift; /* shifted energy */
  return ep + lj->epot_tail; /* unshifted energy */
}

/* energy evaluation, do not change members of `lj' */
INLINE real lj_energyx3d(lj_t *lj, rv3_t *x, real *vir, int *iep, real *ep0, real *eps, real *lap)
{
  real u;
  if (lj->usesq) {
    *iep = lj_energysq3d(lj, x);
    u = (real) (*iep);
  }
  else if (lj->usesw) u = lj_energysw3d(lj, x, vir, lap);
  else u = lj_energylj3d(lj, x, vir, ep0, eps);
  if (lj->usesq || lj->usesw) {
    if (ep0) *ep0 = u;
    if (eps) *eps = u;
  }
  return u;
}

/* compute the energy of the current configuration and set lj->epot */
INLINE real lj_energy3d(lj_t *lj)
{
  return lj->epot = lj_energyx3d(lj, (rv3_t *) lj->x, &lj->vir, &lj->iepot, &lj->epot0, &lj->epots, &lj->lap);
}

real lj_energy(lj_t *lj)
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
        if (lj->usesw & LJ_SWALLPAIRS) { /* save out-of-range pairs, so we can later
            * locate the pair from indices i and j using getpairindex() */
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
      if (laplace)
        lap += 2.f*((672 - 48*3) * dr6 - (192 - 24*3)) * dr6 * dr2;      
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
  if (laplace) *laplace = lap;
  if (f2) for (*f2 = 0.0, i = 0; i < n; i++) *f2 += rv3_sqr(f[i]);  
  return ep + lj->epot_tail; /* unshifted energy */
}

INLINE real lj_force3d(lj_t *lj)
{
  if (lj->usesq) return lj->epot = lj->epot0 = lj->epots = lj_energysq3d(lj,
    (rv3_t *) lj->x); /* no force for square well */
  if (lj->usesw) return lj->epot = lj->epot0 = lj->epots = lj_forcesw3d(lj,
    (rv3_t *) lj->x, (rv3_t *) lj->f,
    lj->pr, &lj->npr, &lj->vir, &lj->f2, &lj->lap);
  return lj->epot = lj_forcelj3d(lj, (rv3_t *) lj->x, (rv3_t *) lj->f,
    &lj->vir, &lj->epot0, &lj->epots, &lj->f2, &lj->lap);
}

real lj_force(lj_t *lj)
{
  return (lj->d == 2) ? lj_force2d(lj) : lj_force3d(lj);
}

/* randomly displace particle i with random amplitude */
INLINE int lj_randmv3d(lj_t *lj, real *xi, real amp)
{
  int i, d;

  i = (int)(rnd0() * lj->n);
  rv3_copy(xi, lj->x + i*3);
  for (d = 0; d < 3; d++) /* displacement */
    xi[d] += (real)(amp * (2.*rnd0() - 1.));
  return i;
}

/* compute energy data for a 3D pair with a square well potential */
INLINE int lj_pairsq3d(real *xi, real *xj, real l, real ra2, real rb2)
{
  real dx[3], dr2;
  dr2 = lj_pbcdist2_3d(dx, xi, xj, l);
  if (dr2 < ra2) return -1000000;
  else if (dr2 < rb2) return 1;
  else return 0;
}

/* return the energy change (square well) from displacing x[i] to xi */
INLINE int lj_depotsq3d(lj_t *lj, int i, real *xi)
{
  int j, n = lj->n, npr = 0;
  real l = lj->l, ra2 = lj->ra2, rb2 = lj->rb2;
  rv3_t *x = (rv3_t *) lj->x;

  for (j = 0; j < n; j++) { /* pair */
    if (j == i) continue;
    npr -= lj_pairsq3d(x[i], x[j], l, ra2, rb2);
    npr += lj_pairsq3d(xi,   x[j], l, ra2, rb2);
  }
  return -npr; /* increased number of pairs == decreased energy */
}

/* commit a particle displacement for a square wel potential */
INLINE void lj_commitsq3d(lj_t *lj, int i, const real *xi, int du)
{
  rv3_copy(lj->x + i*3, xi);
  lj->iepot += du;
  lj->epot += du;
}

/* metropolis for a square well */
INLINE int lj_metrosq3d(lj_t *lj, real amp, real bet)
{
  int i, du;
  real xi[3];

  i = lj_randmv3d(lj, xi, amp);
  du = lj_depotsq3d(lj, i, xi);
  if (metroacc1(du, bet)) {
    lj_commitsq3d(lj, i, xi, du);
    return 1;
  }
  return 0;
}

/* compute energy data for a particle pair, with switched potential  */
INLINE int lj_pairsw3d(lj_t *lj, real *xi, real *xj, real *u, real *vir)
{
  real dx[3], dr2, dr, fscal, psi, ksi;
  dr2 = lj_pbcdist2_3d(dx, xi, xj, lj->l);
  if (dr2 < lj->rc2) {
    dr = (real) sqrt(dr2);
    *u = lj_potsw(lj, dr, &fscal, &psi, &ksi);
    *vir = fscal * dr2; /* f.r */
    return 1;
  } else return 0;
}

/* return the energy change from displacing x[i] to xi */
INLINE real lj_depotsw3d(lj_t *lj, int i, real *xi, real *vir)
{
  int j, n = lj->n;
  real u = 0.f, du = 0.f, dvir = 0.f;
  rv3_t *x = (rv3_t *) lj->x;

  *vir = 0.0f;
  for (j = 0; j < n; j++) { /* pair */
    if (j == i) continue;
    if (lj_pairsw3d(lj, x[i], x[j], &du, &dvir)) {
      u -= du;
      *vir -= dvir;
    }
    if (lj_pairsw3d(lj, xi, x[j], &du, &dvir)) {
      u += du;
      *vir += dvir;
    }
  }
  return u;
}

/* commit a particle displacement 
 * like energysw3d, it does not set pair data, lj->pr
 * call lj_forcesw3d() if it is needed */
INLINE void lj_commitsw3d(lj_t *lj, int i, const real *xi, real du, real dvir)
{
  rv3_copy(lj->x + i*3, xi);
  lj->epot += du;
  lj->vir += dvir;
}

/* Metropolis algorithm */
INLINE int lj_metrosw3d(lj_t *lj, real amp, real bet)
{
  int i;
  real xi[3], du, dvir;

  i = lj_randmv3d(lj, xi, amp);
  du = lj_depotsw3d(lj, i, xi, &dvir);
  if (metroacc1(du, bet)) {
    lj_commitsw3d(lj, i, xi, du, dvir);
    return 1;
  }
  return 0;
}

/* compute energy data for a 3D Lennard-Jones pair */
INLINE int lj_pairlj3d(real *xi, real *xj, real l, real rc2,
    real *u, real *vir)
{
  real dx[3], dr2, invr2, invr6;
  dr2 = lj_pbcdist2_3d(dx, xi, xj, l);
  if (dr2 < rc2) {
    invr2 = 1.0f / dr2;
    invr6 = invr2 * invr2 * invr2;
    *vir = invr6 * (48.f * invr6 - 24.f); /* f.r */
    *u  = 4.f * invr6 * (invr6 - 1.f);
    return 1;
  } else return 0;
}

/* return the energy change from displacing x[i] to xi */
INLINE real lj_depotlj3d(lj_t *lj, int i, real *xi, real *vir)
{
  int j, n = lj->n;
  real l = lj->l, rc2 = lj->rc2, u = 0.f, du = 0.f, dvir = 0.f;
  rv3_t *x = (rv3_t *) lj->x;

  *vir = 0.0f;
  for (j = 0; j < n; j++) { /* pair */
    if (j == i) continue;
    if (lj_pairlj3d(x[i], x[j], l, rc2, &du, &dvir)) {
      u -= du;
      *vir -= dvir;
    }
    if (lj_pairlj3d(xi, x[j], l, rc2, &du, &dvir)) {
      u += du;
      *vir += dvir;
    }
  }
  return u;
}

/* commit a particle displacement */
INLINE void lj_commitlj3d(lj_t *lj, int i, const real *xi, real du, real dvir)
{
  rv3_copy(lj->x + i*3, xi);
  lj->epot += du;
  lj->vir += dvir;
}

INLINE int lj_metrolj3d(lj_t *lj, real amp, real bet)
{
  int i;
  real xi[3], du = 0.f, dvir = 0.f;

  i = lj_randmv3d(lj, xi, amp);
  du = lj_depotlj3d(lj, i, xi, &dvir);
  if (metroacc1(du, bet)) {
    lj_commitlj3d(lj, i, xi, du, dvir);
    return 1;
  }
  return 0;
}

/* return the pair energy between two particles at xi and xj */
INLINE int lj_pair3d(lj_t *lj, real *xi, real *xj, real *u, real *vir)
{
  if (lj->usesq) return lj_pairsq3d(xi, xj, lj->l, lj->ra2, lj->rb2);
  if (lj->usesw) return lj_pairsw3d(lj, xi, xj, u, vir);
  return lj_pairlj3d(xi, xj, lj->l, lj->rc2, u, vir);
}

/* return the energy change from displacing x[i] to xi */
INLINE real lj_depot3d(lj_t *lj, int i, real *xi, real *vir)
{
  if (lj->usesq) return lj_depotsq3d(lj, i, xi);
  if (lj->usesw) return lj_depotsw3d(lj, i, xi, vir);
  return lj_depotlj3d(lj, i, xi, vir);
}

/* this is defined as a macro for du is `int' in sq case, but real in other cases */
#define lj_commit3d(lj, i, xi, du, dvir) { \
  (lj->usesq) ? lj_commitsq3d(lj, i, xi, du) : \
  (lj->usesw) ? lj_commitsw3d(lj, i, xi, du, dvir) : \
  lj_commitlj3d(lj, i, xi, du, dvir); }

/* commit a move */
#define  lj_commit(lj, i, xi, du, dvir) \
  (lj->d == 2 ? lj_commit2d(lj, i, xi, du, dvir) : lj_commit3d(lj, i, xi, du, dvir); )

INLINE int lj_metro3d(lj_t *lj, real amp, real bet)
{
  if (lj->usesq) return lj_metrosq3d(lj, amp, bet);
  if (lj->usesw) return lj_metrosw3d(lj, amp, bet);
  return lj_metrolj3d(lj, amp, bet);
}

/* metropolis algorithm */
INLINE int lj_metro(lj_t *lj, real amp, real bet)
{ return lj->d == 2 ? lj_metro2d(lj, amp, bet) : lj_metro3d(lj, amp, bet); }

/* return the energy change of locally displacing a single atom */
INLINE real lj_dupertl3d(lj_t *lj, real amp)
{
  real dvir, xi[3];
  int i;

  i = lj_randmv3d(lj, xi, amp);
  return lj_depot3d(lj, i, xi, &dvir);
}

/* return the energy change by random displacements of all atoms */
INLINE real lj_dupertg3d(lj_t *lj, real amp)
{
  int i, d, iep;
  rv3_t *nx;
  real du, vir, ep0, eps, lap;

  xnew(nx, lj->n);
  for (i = 0; i < lj->n; i++)
    for (d = 0; d < 3; d++)
      nx[i][d] = lj->x[i*3 + d] + amp * (2*rnd0() - 1);
  du = lj_energyx3d(lj, nx, &vir, &iep, &ep0, &eps, &lap) - lj->epot;
  free(nx);
  return du;
}



/* velocity verlet */
void lj_vv(lj_t *lj, real dt)
{
  int i, nd = lj->n*lj->d;
  real dtl = dt/lj->l;

  for (i = 0; i < nd; i++) { /* VV part 1 */
    lj->v[i] += lj->f[i] * dt *.5f;
    lj->x[i] += lj->v[i] * dtl;
  }
  lj_force(lj); /* calculate the new force */
  for (i = 0; i < nd; i++) /* VV part 2 */
    lj->v[i] += lj->f[i] * dt * .5f;

  lj->ekin = md_ekin(lj->v, nd, lj->dof, &lj->tkin);
  lj->t += dt;
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

  if(udb) for (i = 0; i < n; i++) rv3_zero(h[i]);  /* g * grad g */

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

/* compute volume change, use amp 0.01 ~ 0.02 for 256 system */
INLINE int lj_volmove(lj_t *lj, real amp, real tp, real p)
{
  real lo, ln, loglo, logln, vo, vn, epo, bet = 1.f/tp;
  double dex;

  lo = lj->l;
  vo = lj->vol;
  loglo = (real) log(lo);
  logln = (real) (loglo + amp * (2.f * rnd0() - 1.));
  ln = (real) exp(logln);
  if (ln < lj->rc * 2) return 0; /* box too small */
  epo = lj->epot;
  vn = (real) pow(ln, lj->d);
  lj_setrho(lj, lj->n/vn); /* commit to the new box */
  lj_force(lj);
  dex = bet*(lj->epot - epo + p*(vn - vo)) - (lj->n + 1) * lj->d * (logln - loglo);
/*
  printf("vol %g --> %g, ep %g --> %g, r %g\n", vo, vn, epo, lj->epot, r);
*/
  if (metroacc1(dex, 1.0)) {
    return 1;
  } else {
    lj_setrho(lj, lj->n/vo);
    lj_force(lj); /* inefficient */
    return 0;
  }
}

/* compute volume change, dt < 3e-4 for 256 system
 * cautious, testing version */
INLINE int lj_lgvvolmove(lj_t *lj, real dt, real tp, real p, real dlogvmax)
{
  int d;
  real r, logvo, logvn, dlogv, ln, vn, bet = 1.f/tp;

  logvo = (real) log(lj->vol);
  dlogv = dt * (lj->n + 1 + bet*lj->vir/3.f - bet*p*lj->vol);
  r = (real) grand0();
  dlogv += (real) (r * sqrt(2 * dt));
  dlogv = (real) dblconfine(dlogv, -dlogvmax, dlogvmax); /* change at most 50% */
  logvn = logvo + dlogv;
  ln = (real) exp(logvn/3);
  if (ln < lj->rc * 2) return 0; /* box too small */
  for (vn = 1., d = 0; d < lj->d; d++) vn *= ln;
/*
  printf("pvir %g, bp %g, dlogv %g, rho %g -> %g, p %g\n", (lj->n + 1 + bet*lj->vir/3)/lj->vol, bet*p, dlogv, lj->n/lj->vol, lj->n/vn, lj_calcp(lj, tp));
*/
  lj_setrho(lj, lj->n/vn); /* commit to the new box */
  lj_force(lj);
  /* safety check! */
  return 1;
}


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

  if (lj->d == 3) lj_initfcc3d(lj); else lj_initfcc2d(lj);

  /* init. random velocities */
  for (i = 0; i < n * d; i++) lj->v[i] = (real) (rnd0() - .5);

  lj_shiftcom(lj, lj->v);
  lj_shiftang(lj, lj->x, lj->v);

  lj_force(lj);
  return lj;
}

void lj_close(lj_t *lj)
{
  if (lj->pr) free(lj->pr);
  if (lj->gdg) free(lj->gdg);
  free(lj->x);
  free(lj->v);
  free(lj->f);
  free(lj);
}

/* write position (and velocity)
 * Note 1: *actual* position, not unit position is used
 * Note 2: coordinates are *not* wrapped back into the box */
int lj_writepos(lj_t *lj, const real *x, const real *v, const char *fn)
{
  FILE *fp;
  int i, j, d = lj->d, n = lj->n;
  real l = lj->l;

  if (fn == NULL) fn = "lj.pos";
  xfopen(fp, fn, "w", return -1);

  fprintf(fp, "# %d %d %d %.14e\n", d, lj->n, (v != NULL), l);
  for (i = 0; i < n; i++) {
    for (j = 0; j < d; j++) fprintf(fp, "%16.14f ", x[i*d + j] * l);
    if (v)
      for (j = 0; j < d; j++) fprintf(fp, "%16.14f ", v[i*d + j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* read position file (which may include velocity) */
int lj_readpos(lj_t *lj, real *x, real *v, const char *fn, unsigned flags)
{
  char s[1024], *p;
  FILE *fp;
  int i, j, hasv = 0, next, d = lj->d, n = lj->n;
  const char *fmt;
  real l = lj->l, xtmp;
  double l0;

  if (fn == NULL) fn = "lj.pos";
  xfopen(fp, fn, "r", return -1);

  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') { /* simplified format, old version */
    fprintf(stderr, "Warning: %s has no information line\n", fn);
    rewind(fp);
  } else {
    if (4 != sscanf(s + 1, "%d%d%d%lf", &i, &j, &hasv, &l0)
        || i != d || j != lj->n) {
      fprintf(stderr, "first line is corrupted:\n%s", s);
      goto ERR;
    }
    if (fabs(l0 - l) > 1e-5*l) { /* verify the box size */
      if (flags & LJ_LOADBOX) {
        lj->l = l = (real) l0;
      } else {
        fprintf(stderr, "box mismatch l %g, should be %g\n", l0, l);
        goto ERR;
      }
    }
  }

  fmt = (sizeof(double) == sizeof(real)) ? "%lf%n" : "%f%n";
  for (i = 0; i < n; i++) {
    if (fgets(s, sizeof s, fp) == NULL) goto ERR;
    if (strlen(s) < 10) goto ERR;
    for (p = s, j = 0; j < d; j++, p += next) {
      if (1 != sscanf(p, fmt, &xtmp, &next)) {
        fprintf(stderr, "cannot read i = %d, j = %d\n", i, j);
        goto ERR;
      }
      x[i*d + j] = xtmp / l;
    }
    if (hasv && v != NULL) {
      for (j = 0; j < d; j++, p += next)
        if (1 != sscanf(p, fmt, v + i*d + j, &next)) {
          fprintf(stderr, "cannot read i = %d, j = %d\n", i, j);
          goto ERR;
        }
    }
  }
  fclose(fp);
  return 0;

ERR:
  fprintf(stderr, "position file [%s] appears to be broken on line %d!\n%s\n", fn, i, s);
  fclose(fp);
  return 1;
}

#endif

