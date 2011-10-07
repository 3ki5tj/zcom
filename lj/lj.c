#include "util.c"
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

/* create an open structure */
lj_t *lj_open(int n, int d, real rho, real rcdef, real tp)
{
  lj_t *lj;
  int i;

  xnew(lj, 1);
  lj->n = n;
  lj->d = d;
  lj->dof = n*d - d*(d+1)/2;
  lj->tp = tp;
  xnew(lj->f, n*d);
  xnew(lj->v, n*d);
  xnew(lj->x, n*d);

  lj->rcdef = rcdef;
  lj_setrho(lj, rho);

  if (lj->d == 3) lj_initfcc3d(lj); else lj_initfcc2d(lj);

  /* init. random velocities */
  for (i = 0; i < n*d; i++) lj->v[i] = rnd0() - .5;

  md_shiftcom(lj->x, n, d);
  md_shiftang(lj->x, lj->v, n, d);
  return lj;
}

void lj_close(lj_t *lj)
{
  free(lj->x);
  free(lj->v);
  free(lj->f);
  free(lj);
}

INLINE real lj_pbc(real x, real L)
  { return (real)((x - ((int)((x)+1000.5) - 1000))*L) ; }

INLINE real *lj_vpbc2d(real *v, real L)
  { v[0] = lj_pbc(v[0], L); v[1] = lj_pbc(v[1], L); return v; }

INLINE real lj_pbcdist2_2d(real * RESTRICT dx, const real * RESTRICT a, const real * RESTRICT b, real L)
  { return rv2_sqr(lj_vpbc2d(rv3_diff(dx, a, b), L)); }

INLINE real *lj_vpbc3d(real *v, real L)
  { v[0] = lj_pbc(v[0], L); v[1] = lj_pbc(v[1], L); v[2] = lj_pbc(v[2], L); return v; }

INLINE real lj_pbcdist2_3d(real * RESTRICT dx, const real * RESTRICT a, const real * RESTRICT b, real L)
  { return rv3_sqr(lj_vpbc3d(rv3_diff(dx, a, b), L)); }

static void lj_force2d(lj_t *lj)
{
  real dx[2], fi[2], dr2, dr6, fs, tmp, U, Vir, L = lj->l, rc2 = lj->rc*lj->rc;
  int i, j, d, prcnt = 0, n = lj->n;

  for (i = 0; i < lj->n*2; i++) lj->f[i] = 0;
  for (U = Vir = 0, i = 0; i < n - 1; i++) {
    fi[0] = fi[1] = 0;
    for (j = i+1; j < n; j++) {
      dr2 = lj_pbcdist2_2d(dx, lj->x+i*2, lj->x+j*2, L);
      if (dr2 > rc2) continue;
      dr2 = 1.f/dr2;
      dr6 = dr2*dr2*dr2;
      fs = dr6*(48.f*dr6-24.f); /* f.r */
      Vir += fs;
      fs *= dr2;
      for (d = 0; d < 2; d++) {
        tmp = dx[d]*fs;
        fi[d] += tmp;
        lj->f[j*2+d] -= tmp;
      }
      U += 4*dr6*(dr6-1);
      prcnt++;
    }
    lj->f[i*2+0] += fi[0];
    lj->f[i*2+1] += fi[1];
  }
  lj->epots = U - prcnt*lj->epot_shift; /* shifted energy */
  lj->epot += lj->epot_tail; /* unshifted energy */
  lj->vir = Vir;
  lj->p = lj->rho*lj->tp+lj->vir/(lj->d*lj->vol)+lj->p_tail;
}

/* compute force and virial, return energy */
static void lj_force3d(lj_t *lj)
{
  real dx[3], fi[3], dr2, dr6, fs, tmp, U, Vir, L = lj->l, rc2 = lj->rc*lj->rc;
  int i, j, d, prcnt = 0, n = lj->n;
  real * RESTRICT f = lj->f, * RESTRICT x = lj->x;

  for (i = 0; i < lj->n*3; i++) f[i] = 0;
  for (U = Vir = 0, i = 0; i < n - 1; i++) {
    fi[0] = 0; fi[1] = 0; fi[2] = 0;
    for (j = i+1; j < n; j++) {
      dr2 = lj_pbcdist2_3d(dx, x+i*3, x+j*3, L);
      if (dr2 > rc2) continue;
      dr2 = 1.f/dr2;
      dr6 = dr2*dr2*dr2;
      fs = dr6*(48.f*dr6-24.f); /* f.r */
      Vir += fs;
      fs *= dr2;
      for (d = 0; d < 3; d++) {
        tmp = dx[d]*fs;
        fi[d] += tmp;
        f[j*3+d] -= tmp;
      }
      U += 4*dr6*(dr6-1);
      prcnt++;
    }
    f[i*3+0] += fi[0];
    f[i*3+1] += fi[1];
    f[i*3+2] += fi[2];
  }
  lj->epots = U - prcnt*lj->epot_shift; /* shifted energy */
  lj->epot =  U + lj->epot_tail; /* unshifted energy */
  lj->vir = Vir;
  lj->p = lj->rho*lj->tp+lj->vir/(lj->d*lj->vol)+lj->p_tail;
}

void lj_force(lj_t *lj)
{
  if (lj->d == 3) lj_force3d(lj);
  else if (lj->d == 2) lj_force2d(lj);
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

#endif

