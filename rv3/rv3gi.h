#ifndef RV3GI_H__
#define RV3GI_H__
/* Gauss integral and solid angles */


/* compute the trihedral angle
 * http://planetmath.org/encyclopedia/TrihedralAngle.html
 * tan(omega/2) = vol/den,
 * where
 * den = r1 (r2.r3) + r2 (r1.r3) + r3 (r1.r2) + r1 r2 r3
 */
INLINE real rv3_solidang(const real *v1, const real *v2, const real *v3,
  real * RESTRICT g1, real * RESTRICT g2, real * RESTRICT g3)
{
  real vc[3];
  real r1, r2, r3;
  real ang, vol, den, v2d2, scnum, scden;
  const real eps = (real) 1e-10;

  r1 = rv3_norm(v1);
  r2 = rv3_norm(v2);
  r3 = rv3_norm(v3);
  if (g1 != NULL) rv3_zero(g1);
  if (g2 != NULL) rv3_zero(g2);
  if (g3 != NULL) rv3_zero(g3);

  /* at least two points coincide */
  if (r1 < eps || r2 < eps || r3 < eps) return 0;

  /* the numerator */
  vol = rv3_dot(v3, rv3_cross(vc, v1, v2));

  /* the denominator */
  den  = r1 * rv3_dot(v3, v2);
  den += r2 * rv3_dot(v3, v1);
  den += r3 * (r1*r2 + rv3_dot(v1, v2));

  v2d2 = vol*vol + den*den;

  /* this happens if two vector are opposite
     the solid angle could be evolved from  +/-pi or 0
     but unfortunately, we don't know which one */
  if (v2d2 < eps) return 0;

  if (g1 != NULL && g2 != NULL && g3 != NULL) { /* compute the gradients */
    scnum =  2.f*den/v2d2;
    scden = -2.f*vol/v2d2;

    /* cross products */
    rv3_smul(rv3_cross(g3, v1, v2), scnum);
    rv3_smul(rv3_cross(g1, v2, v3), scnum);
    rv3_smul(rv3_cross(g2, v3, v1), scnum);

    /* compute the contributions to the denominator */
    rv3_lincomb2(vc, v2, v3, r3, r2);
    rv3_sinc(vc, v1, (rv3_dot(v2, v3) + r2*r3)/r1);
    rv3_sinc(g1, vc, scden);

    rv3_lincomb2(vc, v1, v3, r3, r1);
    rv3_sinc(vc, v2, (rv3_dot(v1, v3) + r1*r3)/r2);
    rv3_sinc(g2, vc, scden);

    rv3_lincomb2(vc, v1, v2, r2, r1);
    rv3_sinc(vc, v3, (rv3_dot(v1, v2) + r1*r2)/r3);
    rv3_sinc(g3, vc, scden);
  }

  /* calculate tan(omega/2) */
  ang = (real) atan2(vol, den); /* 0 to pi */
  return 2*ang;
}



/* compute the Gauss integral for the two line segments, with gradients
 *    int_i \int_j (dri X drj).rij/ rij^3,
 * over two line segments rip - ri, rjp - rj, and rij = ri - rj
 * (-2 pi, 2 pi)
 * Note the sign is opposite to that of the dihedral */
INLINE real rv3_solidang2g(const real *ri, const real *rip,
    const real *rj, const real *rjp,
    real * RESTRICT gi, real * RESTRICT gip,
    real * RESTRICT gj, real * RESTRICT gjp)
{
  rv3_t v0, v1, v2, v3, g0, g1, g2, g3, g4, g5;
  real ang1, ang2;

  rv3_diff(v0, ri, rj);
  rv3_diff(v1, ri, rjp);
  rv3_diff(v2, rip, rj);
  rv3_diff(v3, rip, rjp);

  ang1 = rv3_solidang(v0, v1, v2, g0, g1, g2);
  ang2 = rv3_solidang(v2, v1, v3, g3, g4, g5);

  rv3_inc(rv3_inc(rv3_copy(gi,  g0), g1), g4);
  rv3_inc(rv3_inc(rv3_copy(gip, g2), g3), g5);
  rv3_neg(rv3_inc(rv3_inc(rv3_copy(gj,  g0), g2), g3));
  rv3_neg(rv3_inc(rv3_inc(rv3_copy(gjp, g1), g4), g5));

  return ang1 + ang2;
}


/* compute the double integral, old code
 *    \int_i \int_j (dri X drj).rij/ rij^3,
 * over two line segments rip - ri, rjp - rj, and rij = ri - rj
 * (-2 pi, 2 pi)
 * Note the sign is opposite to that of the dihedral */
INLINE real rv3_solidang2(const real *ri, const real *rip,
    const real *rj, const real *rjp)
{
  real v0[3], v1[3], v2[3], v3[3], vc[3];
  double r0, r1, r2, r3;
  double ang, vol, dn1, dn2, dn, tmp;

  r0 = rv3_norm(rv3_diff(v0, ri, rj));
  r1 = rv3_norm(rv3_diff(v1, ri, rjp));
  r2 = rv3_norm(rv3_diff(v2, rip, rj));
  r3 = rv3_norm(rv3_diff(v3, rip, rjp));

  /* avoid coplanar vectors */
  vol = rv3_dot(v0, rv3_cross(vc, v1, v2));
  if (fabs(vol) < 1e-28) return 0;

  /* calculate the denominator */
  tmp = r1*r2 + rv3_dot(v1, v2);
  /* http://planetmath.org/encyclopedia/TrihedralAngle.html
   * tan(omega/2) = vol/den,
   * where
   * den = r1 (r2.r3) + r2 (r1.r3) + r3 (r1.r2) + r1 r2 r3
   * */
  dn1  = r1 * rv3_dot(v0, v2);
  dn1 += r2 * rv3_dot(v0, v1) + r0 * tmp;
  dn2  = r1 * rv3_dot(v3, v2);
  dn2 += r2 * rv3_dot(v3, v1) + r3 * tmp;

  /* calculate tan(omega1/2 + omega2/2) */
  dn = (dn1 + dn2)/(dn1*dn2 - vol*vol);
  ang = atan(fabs(vol) * dn) + (dn < 0 ? M_PI : 0); /* 0 to pi */

  return (real) (vol > 0 ? 2*ang : -2*ang);
}
#endif

