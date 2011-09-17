#include "rng.c"
#ifndef MD_C__
#define MD_C__
#include "md.h"

/* shift the center of mass to zero */
void md_shiftcomw(real *x, real *w, int d, int n)
{
  int i, j;
  real rc, wtot = 0;
  
  if (w) for (i = 0; i < n; i++) wtot += w[i];
  else wtot = n;
  for (j = 0; j < d; j++) {
    for (rc = 0, i = 0; i < n; i++)
      rc += x[i*d+j]*(w ? w[i] : 1.f);
    rc /= wtot;
    for (i = 0; i < n; i++)
      x[i*d+j] -= rc;
  }
}

/* annihilate angular momentum 2d */
void md_shiftang2d(rv2_t *x, rv2_t *v, int n)
{
  int i;
  real am, r2, rp[2];

  for (am = r2 = 0.f, i = 0; i < n; i++) {
    am += rv2_cross(x[i], v[i]);
    r2 += rv2_sqr(x[i]);
  }
  am = -am/r2;
  for (i = 0; i < n; i++) {
    rp[0] = -am*x[i][1];
    rp[1] =  am*x[i][0];
    rv2_inc(v[i], rp);
  }
}

/* annihilate angular momentum 3d
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  / 
 * use
 *    v = c X r
 *   */

void md_shiftang3d(rv3_t *x, rv3_t *v, int n)
{
  int i;
  real ang[3], am[3], dv[3], mat[3][3], inv[3][3];
  const real *xi;
  real xx = 0.f, yy = 0.f, zz = 0.f, xy = 0.f, zx = 0.f, yz = 0.f;

  rv3_zero(am);
  for (i = 0; i < n; i++) {
    rv3_cross(ang, x[i], v[i]);
    rv3_inc(am, ang);
    xi = x[i];
    xx += xi[0]*xi[0];
    yy += xi[1]*xi[1];
    zz += xi[2]*xi[2];
    xy += xi[0]*xi[1];
    yz += xi[1]*xi[2];
    zx += xi[2]*xi[0];
  }
  mat[0][0] = yy+zz;
  mat[1][1] = xx+zz;
  mat[2][2] = xx+yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;
  mat3_inv(inv, mat);
  ang[0] = -rv3_dot(inv[0], am);
  ang[1] = -rv3_dot(inv[1], am);
  ang[2] = -rv3_dot(inv[2], am);
  /* ang is the solution of M^(-1) * I */
  for (i = 0; i < n; i++) {
    rv3_cross(dv, ang, x[i]);
    rv3_inc(v[i], dv);
  }
}

/* return kinetic energy */
real md_ekin(real *v, int nd, int dof, real *tkin)
{
  int i;
  real ekin;

  for (ekin = 0, i = 0; i < nd; i++) ekin += v[i]*v[i];
  if (tkin) *tkin = ekin/dof;
  return ekin *= .5f;
}

/* velocity rescaling thermostat */
void md_vrescale(real tp, real dt, real *ekin, real *tkin, 
    real *v, int nd, int dof)
{
  int i;
  real ekav = .5f*tp*dof, ek1 = *ekin, ek2, s;
  double amp;

  amp = 2*sqrt(ek1*ekav*dt/dof);
  ek2 = ek1 + (ekav - ek1)*dt + (real)(amp*grand0());
  if (ek2 < 0) ek2 = 0;
  s = (real)sqrt(ek2/ek1);
  for (i = 0; i < nd; i++)
    v[i] *= s;
  *ekin = ek2;
  if (tkin) *tkin *= s*s;
}

#endif

