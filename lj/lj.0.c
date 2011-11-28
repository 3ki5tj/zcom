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
  if ((lj->rc = lj->rcdef) > lj->l*.5) lj->rc = lj->l*.5;
  lj->rc2 = lj->rc * lj->rc;
  irc = 1.0/lj->rc;
  irc3 = irc*irc*irc; irc6 = irc3*irc3;
  lj->epot_shift = 4*irc6*(irc6-1);
  if (lj->d == 3) {
    lj->epot_tail = (real)( 8*M_PI*rho*lj->n/9*(irc6 - 3)*irc3 );
    lj->p_tail = (real)( 32*M_PI*rho*rho/9*(irc6 - 1.5)*irc3 );
  } else if (lj->d == 2) {
    lj->epot_tail = (real) (M_PI*rho*lj->n*(.4*irc6 - 1)*irc3*irc);
    lj->p_tail = (real) (M_PI*rho*rho*(1.6*irc6 - 2)*irc3*irc);
  }
}

/* initialize coefficients for the switched potential */
void lj_initsw(lj_t *lj, real rs)
{
  real rs2, rs3, rs6, rs15, dr, dr2, dr3, dr4, f1, f2, f26, f13;

  lj->rs = rs;
  dr = lj->rs - lj->rc;

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
    invr2 = 1.f / (r*r);
    invr6 = invr2 * invr2 * invr2;
    invr8 = invr6 * invr2;
    *fscal = (48.f * invr6 - 24.f) * invr8;
    *psi = (672.f * invr6 - 192.f) * invr8 * invr2;
    *xi = -(10752.f * invr6 - 1920.f) * invr6 * invr6;
    return 4.f * invr6 * (invr6 - 1.f);
  } else if (r < lj->rc) { /* polynomial */
    real dr, dr2, dr3, fs, ps, xs, invr, invr2;
    real a4 = lj->a4, a5 = lj->a5, a6 = lj->a6, a7 = lj->a7;
    invr = 1.f/r;
    dr = r - lj->rc;
    invr2 = invr * invr;
    dr2 = dr * dr;
    dr3 = dr2 * dr;
    fs = dr3*(4.f*a4 + dr*(5.f*a5 + dr*(6.f*a6 + dr*7.f*a7)))*invr;
    *fscal = -fs;
    ps = dr2*(12.f*a4 + dr*(20.f*a5 + dr*(30.f*a6 + dr*42.f*a7)));
    *psi = (ps - fs)*invr2;
    xs = dr*(24.f*a4 + dr*(60.f*a5 + dr*(120.f*a6 + dr*210.0*a7)));
    *xi = (xs*invr - 3.f*(*psi))*invr2;
    return (dr2*dr2)*(a4 + dr*(a5 + dr*(a6 + dr*a7)));
  } else { /* out of range */
    *fscal = 0.f;
    *psi = 0.f;
    *xi = 0.f;
    return 0.f;
  }
}

/* initialize a fcc lattice */
static void lj_initfcc2d(lj_t *lj)
{
  int i, j, id, n1, n = lj->n;
  real a;

  n1 = (int) (pow(2*n, 1.0/lj->d) + .999999); /* # of particles per side */
  a = 1./n1;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++) {
      if ((i+j) % 2 != 0) continue;
      lj->x[id*2 + 0] = (i+.5)*a;
      lj->x[id*2 + 1] = (j+.5)*a;
      id++;
    }
}

/* initialize a fcc lattice */
static void lj_initfcc3d(lj_t *lj)
{
  int i, j, k, id, n1, n = lj->n;
  real a;

  n1 = (int) (pow(2*n, 1.0/lj->d) + .999999); /* # of particles per side */
  a = 1./n1;
  for (id = 0, i = 0; i < n1 && id < n; i++)
    for (j = 0; j < n1 && id < n; j++)
      for (k = 0; k < n1 && id < n; k++) {
        if ((i+j+k) % 2 != 0) continue;
        lj->x[id*3 + 0] = (i+.5)*a;
        lj->x[id*3 + 1] = (j+.5)*a;
        lj->x[id*3 + 2] = (k+.5)*a;
        id++;
      }
}

INLINE real lj_pbc(real x, real L)
  { return (real)((x - ((int)((x)+1000.5) - 1000))*L) ; }

INLINE real *lj_vpbc2d(real *v, real L)
  { v[0] = lj_pbc(v[0], L); v[1] = lj_pbc(v[1], L); return v; }

INLINE real *lj_vpbc3d(real *v, real L)
  { v[0] = lj_pbc(v[0], L); v[1] = lj_pbc(v[1], L); v[2] = lj_pbc(v[2], L); return v; }

INLINE real lj_pbcdist2_3d(real *dx, const real *a, const real *b, real L)
  { return rv3_sqr(lj_vpbc3d(rv3_diff(dx, a, b), L)); }

/* compute energy data for a 3D Lennard-Jones pair */
INLINE int lj_pair3d(real *xi, real *xj, real l, real rc2,
    real *u, real *vir)
{
  real dx[3], dr2, invr2, invr6;
  dr2 = lj_pbcdist2_3d(dx, xi, xj, l);
  if (dr2 < rc2) {
    invr2 = 1.0f / dr2;
    invr6 = invr2 * invr2 * invr2;
    *vir = invr6 * (48.f * invr6 - 24.f);
    *u  = 4.f * invr6 * (invr6 - 1.f);
    return 1;
  } else return 0;
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

/* return the energy change by displacing x[i] to xi */
INLINE real lj_depot3d(lj_t *lj, int i, real *xi, real *vir)
{
  int j, n = lj->n;
  real l = lj->l, rc2 = lj->rc2, u = 0.f, du, dvir;

  *vir = 0.0f;
  for (j = 0; j < n; j++) { /* pair */
    if (j == i) continue;
    if (lj_pair3d(lj->x + i*3, lj->x + j*3, l, rc2, &du, &dvir)) {
      u -= du;
      *vir -= dvir;
    }
    if (lj_pair3d(xi, lj->x + j*3, l, rc2, &du, &dvir)) {
      u += du;
      *vir += dvir;
    }
  }
  return u;
}

/* commit a particle displacement */
INLINE void lj_commit3d(lj_t *lj, int i, const real *xi, real du, real dvir)
{
  rv3_copy(lj->x + i*3, xi);
  lj->epot += du;
  lj->vir += dvir; 
}

INLINE int lj_metro3d(lj_t *lj, real amp, real bet)
{
  int i;
  real xi[3], du, dvir;
  
  i = lj_randmv3d(lj, xi, amp);
  du = lj_depot3d(lj, i, xi, &dvir);
  if (du < 0 || rnd0() < exp(-bet*du)) {
    lj_commit3d(lj, i, xi, du, dvir);
    return 1;
  }
  return 0;
}

/* 3D metropolis algorithm */
INLINE int lj_metro(lj_t *lj, real amp, real bet)
{ return lj->d == 2 ? lj_metro2d(lj, amp, bet) : lj_metro3d(lj, amp, bet); }

/* 3D compute force and virial, return energy */
static void lj_energy3d(lj_t *lj)
{
  real dx[3], dr2, dr6, U, Vir, L = lj->l, rc2 = lj->rc2;
  int i, j, prcnt = 0, n = lj->n;

  for (U = Vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, lj->x + i*3, lj->x + j*3, L);
      if (dr2 > rc2) continue;
      dr2 = 1.f / dr2;
      dr6 = dr2 * dr2 * dr2;
      Vir += dr6 * (48.f * dr6 - 24.f); /* f.r */
      U += 4 * dr6 * (dr6 - 1);
      prcnt++;
    }
  }
  lj->epots = U - prcnt * lj->epot_shift; /* shifted energy */
  lj->epot =  U + lj->epot_tail; /* unshifted energy */
  lj->vir = Vir;
}

void lj_energy(lj_t *lj)
{
  if (lj->d == 2) lj_energy2d(lj); else lj_energy3d(lj);
}

/* 3D compute force and virial, return energy */
static void lj_force3d(lj_t *lj)
{
  real dx[3], fi[3], dr2, dr6, fs, tmp, U, Vir, L = lj->l, rc2 = lj->rc2;
  int i, j, d, prcnt = 0, n = lj->n;
  rv3_t *f = (rv3_t *) lj->f, *x = (rv3_t *) lj->x;

  for (i = 0; i < n; i++) rv3_zero(f[i]);
  for (U = Vir = 0, i = 0; i < n - 1; i++) {
    rv3_zero(fi);
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], L);
      if (dr2 > rc2) continue;
      dr2 = 1.f / dr2;
      dr6 = dr2 * dr2 * dr2;
      fs = dr6 * (48.f * dr6 - 24.f); /* f.r */
      Vir += fs;
      fs *= dr2;
      for (d = 0; d < 3; d++) {
        tmp = dx[d] * fs;
        fi[d] += tmp;
        f[j][d] -= tmp;
      }
      U += 4 * dr6 * (dr6 - 1);
      prcnt++;
    }
    rv3_inc(f[i], fi);
  }
  lj->epots = U - prcnt * lj->epot_shift; /* shifted energy */
  lj->epot =  U + lj->epot_tail; /* unshifted energy */
  lj->vir = Vir;
}

/* compute 3D switched force, save derivative information to lj->pr */
INLINE void lj_forcesw3d(lj_t *lj) 
{
  int i, j, n = lj->n;
  ljpair_t *pr;
  real dx[3], dr2, dr, L = lj->l;
  rv3_t *f = (rv3_t *) lj->f, *x = (rv3_t *) lj->x;
  real fscal, psi, xi;

  lj->epot = lj->lap = lj->vir = 0.f;
  for (i = 0; i < n; i++) rv3_zero(f[i]);
  lj->npr = 0;
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x[i], x[j], L);
      if (dr2 > lj->rc2) continue;

      pr = lj->pr + lj->npr;
      rv3_copy(pr->dx, dx);
      pr->dr2 = dr2;
      dr = (real) sqrt(dr2);
      lj->epot += lj_potsw(lj, dr, &fscal, &psi, &xi);
      pr->phi = -fscal; // phi = u'/r
      pr->psi = psi;  // psi = phi'/r
      pr->xi = xi;  // xi = psi'/r
      lj->lap += 2.f*psi*dr2 - 6.f*fscal;
      lj->vir -= fscal*dr2;
      rv3_smul(dx, fscal);
      pr->i = i;
      pr->j = j;
      rv3_inc(f[i], dx);
      rv3_dec(f[j], dx);
      lj->npr++;
    }
  }
  for (lj->f2 = 0.0, i = 0; i < n; i++) lj->f2 += rv3_sqr(f[i]);
}

void lj_force(lj_t *lj)
{
  if (lj->d == 2) lj_force2d(lj);
  else if (lj->usesw) lj_forcesw3d(lj);
  else lj_force3d(lj);
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
  for (i = 0; i < n * d; i++) lj->v[i] = rnd0() - .5;

  md_shiftcom(lj->x, n, d);
  md_shiftang(lj->x, lj->v, n, d);

  lj_force(lj);
  return lj;
}

void lj_close(lj_t *lj)
{
  free(lj->x);
  free(lj->v);
  free(lj->f);
  free(lj);
}
#endif

