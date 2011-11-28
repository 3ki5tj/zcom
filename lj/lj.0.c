#include "util.h"
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
INLINE int lj_calcpair3d(ljpair_t *pr, real *xi, real *xj, real l, real rc2)
{
  real dx[3], dr2, invr2, invr6;
  dr2 = lj_pbcdist2_3d(dx, xi, xj, l);
  if (dr2 < rc2) {
    invr2 = 1.0f / dr2;
    invr6 = invr2 * invr2 * invr2;
    pr->vir = invr6 * (48.f * invr6 - 24.f);
    pr->u  = 4.f * invr6 * (invr6 - 1.f);
    return pr->in = 1;
  } else return pr->in = 0;
}

INLINE void lj_initmc3d(lj_t *lj)
{
  int i, j, n = lj->n;
  real l = lj->l, rc2 = lj->rc * lj->rc;
  ljpair_t *pr;

  lj->epot = 0.f;
  lj->vir = 0.f;
  for (i = 0; i < n; i++) 
    for (j = i + 1; j < n; j++) {
      pr = lj->pair + i * n + j;
      lj_calcpair3d(pr, lj->x + i*3, lj->x + j*3, l, rc2);
      lj->epot += pr->u;
      lj->vir += pr->vir;
    }
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
  real l = lj->l, rc2 = lj->rc * lj->rc, u = 0.f;
  ljpair_t *opr, *npr;

  *vir = 0.0f;
  for (j = 0; j < n; j++) { /* pair */
    if (j == i) continue;
    opr = (j > i) ? (lj->pair + i * n + j) : (lj->pair + j * n + i);
    if (opr->in) {
      *vir -= opr->vir;
      u -= opr->u;
    }
    npr = lj->npair + j;
    if (lj_calcpair3d(npr, xi, lj->x + j*3, l, rc2)) {
      *vir += npr->vir;
      u += npr->u;
    }
  }
  return u;
}

/* commit a particle displacement */
INLINE void lj_commit3d(lj_t *lj, int i, const real *xi, real du, real dvir)
{
  int j, n = lj->n;
  ljpair_t *opr, *npr;

  rv3_copy(lj->x + i*3, xi);
  for (j = 0; j < n; j++) {
    if (j == i) continue;
    opr = (j > i) ? (lj->pair + i * n + j) : (lj->pair + j * n + i);
    npr = lj->npair + j;
    if ( (opr->in = npr->in) ) {
      opr->u = npr->u;
      opr->vir = npr->vir;
    }
  }
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
  real dx[3], dr2, dr6, U, Vir, L = lj->l, rc2 = lj->rc * lj->rc;
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
  real dx[3], fi[3], dr2, dr6, fs, tmp, U, Vir, L = lj->l, rc2 = lj->rc * lj->rc;
  int i, j, d, prcnt = 0, n = lj->n;
  real * RESTRICT f = lj->f, * RESTRICT x = lj->x;

  for (i = 0; i < n * 3; i++) f[i] = 0.f;
  for (U = Vir = 0, i = 0; i < n - 1; i++) {
    for (d = 0; d < 3; d++) fi[d] = 0.f;
    for (j = i + 1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x + i*3, x + j*3, L);
      if (dr2 > rc2) continue;
      dr2 = 1.f / dr2;
      dr6 = dr2 * dr2 * dr2;
      fs = dr6 * (48.f * dr6 - 24.f); /* f.r */
      Vir += fs;
      fs *= dr2;
      for (d = 0; d < 3; d++) {
        tmp = dx[d] * fs;
        fi[d] += tmp;
        f[j * 3 + d] -= tmp;
      }
      U += 4 * dr6 * (dr6 - 1);
      prcnt++;
    }
    for (d = 0; d < 3; d++) f[i*3 + d] = fi[d];
  }
  lj->epots = U - prcnt * lj->epot_shift; /* shifted energy */
  lj->epot =  U + lj->epot_tail; /* unshifted energy */
  lj->vir = Vir;
}

void lj_force(lj_t *lj)
{
  if (lj->d == 2) lj_force2d(lj); else lj_force3d(lj);
}

/* velocity verlet */
void lj_vv(lj_t *lj, real dt)
{
  int i, nd = lj->n*lj->d;
  real dth = dt*.5f, dtl = dt/lj->l;

  for (i = 0; i < nd; i++) { /* VV part 1 */
    lj->v[i] += lj->f[i]*dth;
    lj->x[i] += lj->v[i]*dtl;
  }
  lj_force(lj); /* calculate the new force */
  for (i = 0; i < nd; i++) /* VV part 2 */
    lj->v[i] += lj->f[i]*dth;
  
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

  xnew(lj->pair, n * n);
  xnew(lj->npair, n);
  if (d == 2) lj_initmc2d(lj); else lj_initmc3d(lj);
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

