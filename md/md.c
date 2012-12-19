#include "rng.c"
#ifndef MD_C__
#define MD_C__
#include "md.h"

/* shift the center of mass to zero */
INLINE void md_shiftcomw(real * RESTRICT x, const real * RESTRICT w, int n, int d)
{
  int i, j;
  real rc, wtot = 0;
  
  if (w) for (i = 0; i < n; i++) wtot += w[i];
  else wtot = (real) n;
  for (j = 0; j < d; j++) {
    for (rc = 0, i = 0; i < n; i++)
      rc += x[i*d+j]*(w ? w[i] : 1.f);
    rc /= wtot;
    for (i = 0; i < n; i++)
      x[i*d+j] -= rc;
  }
}

/* annihilate angular momentum 2d */
INLINE void md_shiftang2d(rv2_t * RESTRICT x, rv2_t * RESTRICT v, int n)
{
  int i;
  real am, r2, xc[2] = {0,0}, xi[2];

  for (i = 0; i < n; i++) rv2_inc(xc, x[i]);
  rv2_smul(xc, 1.f/n);
  for (am = r2 = 0.f, i = 0; i < n; i++) {
    rv2_diff(xi, x[i], xc);
    am += rv2_cross(xi, v[i]);
    r2 += rv2_sqr(x[i]);
  }
  am = -am/r2;
  for (i = 0; i < n; i++) {
    rv2_diff(xi, x[i], xc);
    v[i][0] += -am*xi[1];
    v[i][1] +=  am*xi[0];
  }
}

/* annihilate angular momentum 3d
 * solve
 *   /  y^2 + z^2    -x y      -x y      \
 *   |  -x y       X^2 + z^2   -y z      |  c  =  I
 *   \  -x z         -y z     x^2 + y^2  / 
 * use a velocity field 
 *    v = c X r
 *   */

INLINE void md_shiftang3d(rv3_t *x, rv3_t *v, int n)
{
  int i;
  real xc[3] = {0,0,0}, xi[3], ang[3], am[3] = {0,0,0}, dv[3], mat[3][3], inv[3][3];
  real xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;

  for (i = 0; i < n; i++) rv3_inc(xc, x[i]);
  rv3_smul(xc, 1.f/n);
  for (i = 0; i < n; i++) {
    rv3_diff(xi, x[i], xc);
    rv3_cross(ang, xi, v[i]);
    rv3_inc(am, ang);
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
  rm3_inv(inv, mat);
  ang[0] = -rv3_dot(inv[0], am);
  ang[1] = -rv3_dot(inv[1], am);
  ang[2] = -rv3_dot(inv[2], am);
  /* ang is the solution of M^(-1) * I */
  for (i = 0; i < n; i++) {
    rv3_diff(xi, x[i], xc);
    rv3_cross(dv, ang, xi);
    rv3_inc(v[i], dv);
  }
}

/* return kinetic energy */
INLINE real md_ekin(const real *v, int nd, int dof, real *tkin)
{
  int i;
  real ekin;

  for (ekin = 0, i = 0; i < nd; i++) ekin += v[i]*v[i];
  if (tkin) *tkin = ekin/dof;
  return ekin *= .5f;
}

/* velocity scaling: for regular MD
 * ekt is the time-averaged ek, may not be *ekin  */
INLINE void md_vscale(real *v, int nd, int dof, real tp, real ekt, real *ekin, real *tkin)
{
  int i;
  real ekav = .5f*tp*dof, s;

  s = (real) sqrt( ekav / ekt );
  for (i = 0; i < nd; i++)
    v[i] *= s;
  if (ekin) *ekin *= s*s;
  if (tkin) *tkin *= s*s;
}

/* velocity rescaling thermostat */
INLINE void md_vrescale(real *v, int nd, int dof, real tp, real dt, real *ekin, real *tkin)
{
  int i;
  real ekav = .5f*tp*dof, ek1 = *ekin, ek2, s;
  double amp;

  amp = 2*sqrt(ek1*ekav*dt/dof);
  ek2 = ek1 + (ekav - ek1)*dt + (real)(amp*grand0());
  if (ek2 < 0) ek2 = 0;
  s = (real) sqrt(ek2/ek1);
  for (i = 0; i < nd; i++)
    v[i] *= s;
  *ekin = ek2;
  if (tkin) *tkin *= s*s;
}


/* Exact velocity rescaling thermostat */
INLINE void md_vrescalex(real *v, int nd, int dof, real tp, real dt, real *ekin, real *tkin)
{
  int i;
  real ekav = .5f*tp*dof, ek1 = *ekin, ek2, s, c = 0., r;
  
  if (dt < 10.0) c = exp(-dt);
  r = (real) grand0();
  ek2 = (real)( ek1 + (1.f-c)*(ekav*(randgausssum(dof-1) + r*r)/dof - ek1)
    + 2.0f*r*sqrt(c*(1.f-c)*ekav/dof*ek1) );
  if (ek2 < 0) ek2 = 0;
  s = (real) sqrt(ek2/ek1);
  for (i = 0; i < nd; i++)
    v[i] *= s;
  *ekin = ek2;
  if (tkin) *tkin *= s*s;
}

/* backup thermostat: velocity-rescaling according to a Monte-Carlo move */
INLINE int md_mcvrescale(real *v, int nd, int dof, real tp, real dt, real *ekin, real *tkin)
{
  int i;
  real ek1 = *ekin, ek2, s;
  double logek1, logek2, r;

  logek1 = log(ek1);
  logek2 = logek1 + dt*(2.f*rnd0() - 1);
  ek2 = exp(logek2);
  r = (ek2-ek1)/tp - .5*dof*(logek2 - logek1);
  if (r <= 0 || rnd0() < exp(-r)) {
    s = (real) sqrt(ek2/ek1);
    for (i = 0; i < nd; i++)
      v[i] *= s;
    *ekin = ek2;
    if (tkin) *tkin *= s*s;
    return 1;
  } else { /* do nothing otherwise */
    return 0;
  }
}

INLINE int md_mcvrescale2d(rv2_t * RESTRICT v, int n, int dof, real tp, real dt, real * RESTRICT ekin, real * RESTRICT tkin)
    { return md_mcvrescale((real *) v, n*2, dof, tp, dt, ekin, tkin); }
INLINE int md_mcvrescale3d(rv3_t * RESTRICT v, int n, int dof, real tp, real dt, real * RESTRICT ekin, real * RESTRICT tkin)
    { return md_mcvrescale((real *) v, n*3, dof, tp, dt, ekin, tkin); }


/* Nose-Hoover thermostat */
INLINE void md_hoover(real *v, int nd, int dof, real tp, real dt,
   real *zeta, real Q, real *ekin, real *tkin)
{
  int i;
  real ek1 = *ekin, ek2, s, dt2 = .5f*dt;
  
  *zeta += (2.f*ek1 - dof * tp)/Q*dt2;
  
  s = (real) exp(-(*zeta)*dt);
  for (i = 0; i < nd; i++) v[i] *= s;
  ek2 = ek1 * (s*s);
  *ekin = ek2;
  if (tkin) *tkin *= s*s;
  
  *zeta += (2.f*ek2 - dof * tp)/Q*dt2;
}

INLINE void md_hoover2d(rv2_t *v, int n, int dof, real tp, real dt,
   real *zeta, real Q, real *ekin, real *tkin)
  { md_hoover((real *)v, n*2, dof, tp, dt, zeta, Q, ekin, tkin); }

INLINE void md_hoover3d(rv3_t *v, int n, int dof, real tp, real dt,
   real *zeta, real Q, real *ekin, real *tkin)
  { md_hoover((real *)v, n*3, dof, tp, dt, zeta, Q, ekin, tkin); }

/* Nose-Hoover chain thermostat */
INLINE void md_nhchain(real *v, int nd, int dof, real tp, real scl, real dt,
   real *zeta, const real *Q, int M, real *ekin, real *tkin)
{
  int i, j;
  real ek1 = *ekin, ek2, s, dt2 = .5f*dt, dt4 = .25f*dt, G, xp = 1.f;
 
  /* propagate the chain */
  for (j = M-1; j > 0; j--) {
    if (j < M-1) {
      xp = (real) exp(-dt4*zeta[j+1]);
      zeta[j] *= xp;
    }
    G = (Q[j-1]*zeta[j-1]*zeta[j-1] - tp)/Q[j];
    zeta[j] += G * dt2;
    if (j < M-1)
      zeta[j] *= xp;
  }

  /* the first thermostat variable */
  if (M >= 2) {
    xp = exp(-dt4*zeta[1]);
    zeta[0] *= xp;
  }
  G = (scl * 2.f * ek1 - dof * tp)/Q[0];
  zeta[0] += G * dt2;
  if (M >= 2)
    zeta[0] *= xp;
 
  /* scale the velocities */ 
  s = (real) exp(-(*zeta)*dt);
  for (i = 0; i < nd; i++) v[i] *= s;
  ek2 = ek1 * (s*s);
  *ekin = ek2;
  if (tkin) *tkin *= s*s;
 
  /* the first thermotat variable */ 
  if (M >= 2) {
    xp = exp(-dt4*zeta[1]);
    zeta[0] *= xp;
  }
  G = (scl * 2.f * ek1 - dof * tp)/Q[0];
  zeta[0] += G * dt2;
  if (M >= 2)
    zeta[0] *= xp;
  
  /* propagate the chain */
  for (j = M-1; j > 0; j--) {
    if (j < M-1) {
      xp = (real) exp(-dt4*zeta[j+1]);
      zeta[j] *= xp;
    }
    G = (Q[j-1]*zeta[j-1]*zeta[j-1] - tp)/Q[j];
    zeta[j] += G * dt2;
    if (j < M-1)
      zeta[j] *= xp;
  }
}

INLINE void md_nhchain2d(rv3_t *v, int n, int dof, real tp, real scl, real dt,
   real *zeta, const real *Q, int M, real *ekin, real *tkin)
  { md_nhchain((real *)v, n*2, dof, tp, scl, dt, zeta, Q, M, ekin, tkin); }

INLINE void md_nhchain3d(rv3_t *v, int n, int dof, real tp, real scl, real dt,
   real *zeta, const real *Q, int M, real *ekin, real *tkin)
  { md_nhchain((real *)v, n*3, dof, tp, scl, dt, zeta, Q, M, ekin, tkin); }

/* Anderson thermostat */
INLINE void md_andersen(real *v, int n, int d, real tp)
{
  int i, j;

  tp = sqrt(tp);
  i = (int)(rnd0() * n);
  for (j = 0; j < d; j++) v[i*d + j] = tp * grand0();
}

/* Langevin thermostat */
INLINE void md_langevin(real *v, int n, int d, real tp, real dt)
{
  int i;
  real c = (real) exp(-dt), amp;

  amp = sqrt((1 - c*c) * tp);
  for (i = 0; i < n*d; i++)
    v[i] = c*v[i] + amp*grand0();
}

/* Nose-Hoover thermostat/barostat
 * set cutoff to half of the box */
INLINE void md_hoovertp(real *v, int n, int d, int dof, real dt, 
    real tp, real pext, real *zeta, real *eta, real Q, real W,
    real vol, real vir, real ptail, int ensx,
    real *ekin, real *tkin)
{
  int i;
  real xp, pint, s, dt2 = dt*.5f, dt4 = dt*.25f;

  /* thermostat */
  *zeta += (2.f * (*ekin) + W * (*eta) * (*eta) - (dof + 1) * tp) * dt2/Q;
  xp = (real) exp(-(*zeta)*dt4); /* zeta won't change until the end */

  /* barostat */
  *eta *= xp;
  pint = (vir + 2.f * (*ekin))/ (d * vol) + ptail;
  *eta += ((pint - pext)*vol + (1 - ensx) * tp)*d*dt2/W;
  *eta *= xp;

  /* scaling velocity */
  s = exp( -dt * (*zeta + *eta) );
  for (i = 0; i < d * n; i++) v[i] *= s;
  *ekin *= s*s;
  *tkin *= s*s;

  /* barostat */
  *eta *= xp;
  pint = (vir + 2.f * (*ekin))/ (d * vol) + ptail;
  *eta += ((pint - pext)*vol + (1 - ensx) * tp)*d*dt2/W;
  *eta *= xp;

  /* thermostat */
  *zeta += (2.f * (*ekin) + W * (*eta) * (*eta) - (dof + 1) * tp) * dt2/Q;
}

/* Nose-Hoover chain thermostat/barostat
 * set cutoff to half of the box */
INLINE void md_nhchaintp(real *v, int n, int d, int dof, real dt, 
    real tp, real pext, real *zeta, real *eta, const real *Q, int M, real W,
    real vol, real vir, real ptail, int ensx,
    real *ekin, real *tkin)
{
  int i, j;
  real xpz, pint, s, dt2 = dt*.5f, dt4 = dt*.25f, G, xp;

  /* 1. thermostat */
  /* 1.A propagate the chain */
  for (j = M-1; j > 0; j--) {
    if (j < M-1) {
      xp = (real) exp(-dt4*zeta[j+1]);
      zeta[j] *= xp;
    }
    G = (Q[j-1]*zeta[j-1]*zeta[j-1] - tp)/Q[j];
    zeta[j] += G * dt2;
    if (j < M-1)
      zeta[j] *= xp;
  }

  /* 1.B the first thermostat variable */
  if (M >= 2) {
    xp = exp(-dt4*zeta[1]);
    zeta[0] *= xp;
  }
  G = (2.f * (*ekin) + W * (*eta) * (*eta) - (dof + 1) * tp) / Q[0];
  zeta[0] += G * dt2;
  if (M >= 2)
    zeta[0] *= xp;
  xpz = (real) exp(-zeta[0]*dt4); /* zeta won't change until the end */

  /* 2. barostat */
  *eta *= xpz;
  pint = (vir + 2.f * (*ekin))/ (d * vol) + ptail;
  *eta += ((pint - pext)*vol + (1 - ensx) * tp)*d*dt2/W;
  *eta *= xpz;

  /* 3. scaling velocity */
  s = exp( -dt * (zeta[0] + *eta) );
  for (i = 0; i < d * n; i++) v[i] *= s;
  *ekin *= s*s;
  *tkin *= s*s;

  /* 4. barostat */
  *eta *= xpz;
  pint = (vir + 2.f * (*ekin))/ (d * vol) + ptail;
  *eta += ((pint - pext)*vol + (1 - ensx) * tp)*d*dt2/W;
  *eta *= xpz;

  /* 5. thermostat */
  /* 5.A the first thermotat variable */
  if (M >= 2) {
    xp = exp(-dt4*zeta[1]);
    zeta[0] *= xp;
  }
  G = (2.f * (*ekin) + W * (*eta) * (*eta) - (dof + 1) * tp) / Q[0];
  zeta[0] += G * dt2;
  if (M >= 2)
    zeta[0] *= xp;

  /* 5.B propagate the chain */
  for (j = M-1; j > 0; j--) {
    if (j < M-1) {
      xp = (real) exp(-dt4*zeta[j+1]);
      zeta[j] *= xp;
    }
    G = (Q[j-1]*zeta[j-1]*zeta[j-1] - tp)/Q[j];
    zeta[j] += G * dt2;
    if (j < M-1)
      zeta[j] *= xp;
  }
}

/* Langevin barostat
 *   d eta / dt = -zeta * eta
 *      + [ (Pint - Pext) * V + (1 - ensx) T ] * d / W 
 *      + sqrt( 2 * zeta * T / W ) xi
 * the ideal-gas part of the pressure, Pint, is computed as \sum p^2/m / V
 * the additional volume distribution weight is V^(-ensx)
 * the scaling is r = r*s, p = p/s;
 * set the cutoff rc to half of the box */
INLINE void md_langtp(real *v, int n, int d, real dt, 
    real tp, real pext, real zeta, real *eta, real W,
    real vol, real vir, real ptail, int ensx,
    real *ekin, real *tkin)
{
  int i;
  real xp, pint, s, dt2 = dt*.5f, dt4 = dt*.25f, amp;

  xp = (real) exp(-zeta*dt4);
  amp = (real) sqrt(2.f*zeta*tp/W*dt2);

  /* barostat: first half to update *eta */
  *eta *= xp;
  pint = (vir + 2.f * (*ekin))/ (d * vol) + ptail;
  *eta += ((pint - pext)*vol + (1 - ensx) * tp)*d*dt2/W;
  *eta += amp*grand0(); /* random noise */
  *eta *= xp;

  /* scaling velocity */
  s = exp( -dt * (*eta) );
  for (i = 0; i < d * n; i++) v[i] *= s;
  *ekin *= s*s;
  *tkin *= s*s;

  /* barostat: second half to update *eta */
  *eta *= xp;
  pint = (vir + 2.f * (*ekin))/ (d * vol) + ptail;
  *eta += ((pint - pext)*vol + (1 - ensx) * tp)*d*dt2/W;
  *eta += amp*grand0(); /* random noise */
  *eta *= xp;
}

/* position Langevin barostat,
 * limiting case, zeta -> inf., of the full Langevin barostat
 * barodt = dt/(W*zeta), and d * d(eta) = d lnv
 * the ideal-gas part of the pressure is computed as \sum p^2/m / V
 * the scaling is r = r*s, p = p/s;
 * set cutoff to half of the box */
INLINE void md_langtp0(real *v, int n, int d, real barodt,
   real tp, real pext, real *vol, real vir, real ptail, int ensx,
   real *ekin, real *tkin)
{
  int i;
  real pint, amp, vn, s, dlnv;

  /* compute the internal pressure */
  /* note only with half-box cutoff, the formula is accurate */
  pint = (vir + 2.f * (*ekin))/ (d * (*vol)) + ptail;

  amp = (real) sqrt(2.f * barodt);
  dlnv = ((pint - pext) * (*vol)/tp + 1 - ensx)*barodt + amp*grand0();
  vn = log(*vol) + dlnv;
  vn = exp(vn);

  s = (real) pow(*vol/vn, 1.0/3);
  for (i = 0; i < d * n; i++) v[i] *= s;
  *ekin *= s*s;
  *tkin *= s*s;

  *vol = vn;
}


/* sinc(x) = (e^x - e^(-x))/(2 x) */
INLINE double md_mysinc(double x)
{
  double x2 = x*x;

  if (fabs(x) < 1e-2) /* series expansion */
    return 1 + (1 + (1 + x2/42.0)*x2/20.0)*x2/6.0;
  else
    return .5 * (exp(x) - exp(-x))/x;
}

/* Nose-Hoover position update */
INLINE void md_hoovertpdr(real *r, const real *v, int nd,
    real *xp, real l, real eta, real dt)
{
  int i;
  real dtxp, xph, etadt2;
  
  /* r' = r*exp(eta*dt) + v*(exp(eta*dt) - 1)/eta 
   * now exp(eta*dt) is the volume scaling factor
   * so for the reduced coordinates R = r*exp(-eta*dt)
   * R' = R + v*(1 - exp(-eta*dt))/eta; */
  etadt2 = eta * dt * .5f;
  xph = (real) exp(etadt2);
  *xp = xph * xph;
  dtxp = 1.f/xph * dt * md_mysinc(etadt2) / l;
/*  
  dtxp = (1 - 1/(*xp))/eta/l;
*/  
  for (i = 0; i < nd; i++)
    r[i] += v[i] * dtxp;
}

#endif

