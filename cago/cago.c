#include "util.h"
#include "av.h"
#include "pdb.c"
#include "rng.c"
#include "rotfit.c"
#include "md.c"
#ifndef CAGO_C__
#define CAGO_C__
#include "cago.h"

/* initialize data for the potential energy function
 * compute the reference bonds and angles */
INLINE int cago_initpot(cago_t *go)
{
  int i, j, n = go->n;
  real dr2;

  /* calculate reference geometric parameters */
  xnew(go->bref, n - 1); /* bonds */
  for (i = 0; i < n - 1; i++)
    go->bref[i] = rv3_dist(go->xref[i], go->xref[i+1]);

  xnew(go->aref, n - 2); /* angles */
  for (i = 1; i < n - 1; i++) {
    go->aref[i-1]  = rv3_ang(go->xref[i-1], go->xref[i], go->xref[i+1],
      NULL, NULL, NULL);
  }

  xnew(go->dref, n - 3); /* dihedrals */
  for (i = 0; i < n - 3; i++) {
    go->dref[i] = rv3_calcdih(NULL,
      go->xref[i], go->xref[i+1], go->xref[i+2], go->xref[i+3], 0);
  }

  /* reference pair distances */
  xnew(go->r2ref, n*n);
  for (i = 0; i < n - 1; i++) {
    for (j = i+1; j < n; j++) {
      dr2 = rv3_dist2(go->xref[i], go->xref[j]);
      go->r2ref[j*n + i] = go->r2ref[i*n + j] = dr2;
    }
  }
  return 0;
}

/* return cago_t from pdb file fnpdb
 * `rcc' is the cutoff radius for defining contacts
 * the last fours bits of `flags' is reserved for defining PDB_CONTACT_XXX */
cago_t *cago_openx(const char *fnpdb, real kb, real ka, real kd1, real kd3,
    real nbe, real nbc, real rcc, unsigned flags)
{
  cago_t *go;
  int i, j;
  pdbmodel_t *pm;
  pdbaac_t *c;

  xnew(go, 1);

  /* read all atoms from the .pdb file */
  if ((pm = pdbm_read(fnpdb, 0)) == NULL)
    return NULL;
  go->iscont = pdbm_contact(pm, rcc,
      flags & 0xf, /* the last four bits is translated to PDB_CONTACT_XXX */
      3, /* `a' and `d' in -a-b-c-d- are excluded from being contacts */
      flags & CAGO_VERBOSE);
  /* parse it into amino-acid residues */
  if ((c = pdbaac_parse(pm, 0)) == NULL)
    return NULL;
  pdbm_free(pm);

  go->n = c->nres;
  go->dof = go->n*3 - 6;
  go->kb = kb;
  go->ka = ka;
  go->kd1 = kd1;
  go->kd3 = kd3;
  go->nbe = nbe;
  go->nbc = nbc;
  xnew(go->xref, go->n);
  xnew(go->aa, go->n);
  /* extract the C-alpha coordinates */
  for (i = 0; i < go->n; i++) {
    for (j = 0; j < 3; j++)
      go->xref[i][j] = (real) c->res[i].xca[j];
    go->aa[i] = c->res[i].aa;
  }
  pdbaac_free(c); /* throw away pdbaac */

  /* compute the reference bond length, angles, etc. */
  cago_initpot(go);

  for (go->ncont = 0, i = 0; i < go->n - 1; i++)
    for (j = i+1; j < go->n; j++)
      go->ncont += go->iscont[ i*go->n + j ];
  return go;
}

/* destroy cago_t */
void cago_close(cago_t *go)
{
  if (go->x) free(go->x);
  if (go->v) free(go->v);
  if (go->f) free(go->f);
  if (go->x1) free(go->x1);
  if (go->iscont) free(go->iscont);
  free(go->bref);
  free(go->aref);
  free(go->dref);
  free(go->r2ref);
  free(go->xref);
  free(go->aa);
  memset(go, '\0', sizeof(*go));
  free(go);
}

/* remove center of mass motion, linear and angular */
INLINE void cago_rmcom(cago_t *go, rv3_t *x, rv3_t *v)
{
  md_shiftcom3d(x, go->n);
  md_shiftcom3d(v, go->n);
  md_shiftang3d(x, v, go->n);
}

/* initialize molecular dynamics
 *  o create an initial structure
 *    if rndamp >= 0, start from the reference structure,
 *      with a random disturbance of rndamp
 *    if rndamp < 0, start from a nearly-straight chain, 
 *      with a disturbance of rndamp in the x, y directions
 *  o initialize the velocity with the center of mass motion removed
 *  o compute the initial force and energy
 * */
INLINE int cago_initmd(cago_t *go, double rndamp, double T0)
{
  int i, j, n = go->n;
  real s, dx[3];

  xnew(go->f, n);
  xnew(go->v, n);
  xnew(go->x, n);
  xnew(go->x1, n);

  /* initialize position */
  if (rndamp < 0) { /* open chain */
    rndamp *= -1;
    for (i = 0; i < n-1; i++) {
      for (j = 0; j < 3; j++)
        dx[j] = (j == 0) ? 1.f : rndamp*(2.f*rnd0() - 1);
      rv3_normalize(dx);
      rv3_smul(dx, go->bref[i]);
      rv3_add(go->x[i+1], go->x[i], dx);
    }
  } else { /* copy from xref, slightly disturb it */
    for (i = 0; i < n; i++) {
      rv3_copy(go->x[i], go->xref[i]);
      for (j = 0; j < 3; j++)
        go->x[i][j] += rndamp*(2.f*rnd0() - 1);
    }
  }
  go->epotref = cago_force(go, go->xref, go->f);
  go->epot = cago_force(go, go->x, go->f);

  /* initialize velocities */
  for (j = 0; j < 3; j++)
    for (i = 0; i < n; i++)
      go->v[i][j] = rnd0() - .5;
  cago_rmcom(go, go->x, go->v); /* remove center of mass motion */
  for (s = 0, i = 0; i < n; i++)
    s += rv3_sqr(go->v[i]);
  s = sqrt( (3*n*T0)/s );
  for (i = 0; i < n; i++) {
    rv3_smul(go->v[i], s);
  }
  go->ekin = cago_ekin(go, go->v);
  go->rmsd = cago_rmsd(go, go->x, NULL);
  go->t = 0;
  return 0;
}

/* bond energy 1/2 k (r - r0)^2 */
INLINE real potbond(rv3_t a, rv3_t b, real r0, real k, 
    rv3_t fa, rv3_t fb)
{
  real dx[3], r, dr, amp;
  
  r = rv3_norm( rv3_diff(dx, a, b) );
  dr = r - r0;
  if (fa != NULL) {
    amp = k * dr / r;
    rv3_sinc(fa, dx, -amp);
    rv3_sinc(fb, dx,  amp);
  }
  return .5f * k * dr * dr;
}

/* harmonic angle 1/2 k (ang - ang0)^2 */
INLINE real potang(rv3_t a, rv3_t b, rv3_t c, real ang0, real k,
    rv3_t fa, rv3_t fb, rv3_t fc)
{
  real dang, amp, ga[3], gb[3], gc[3];

  if (fa) { /* compute gradient */ 
    dang = rv3_ang(a, b, c, ga, gb, gc) - ang0;
    amp = -k * dang;
    rv3_sinc(fa, ga, amp);
    rv3_sinc(fb, gb, amp);
    rv3_sinc(fc, gc, amp);
  } else {
    dang = rv3_ang(a, b, c, NULL, NULL, NULL) - ang0;
  }
  return .5f * k * dang * dang;
}

/* 1-3 dihedral: k1 * (1 - cos(dang)) + k3 * (1 - cos(3*dang)) */
INLINE real potdih13(rv3_t a, rv3_t b, rv3_t c, rv3_t d, real ang0,
    real k1, real k3, rv3_t fa, rv3_t fb, rv3_t fc, rv3_t fd) 
{
  real dang, amp, ga[3], gb[3], gc[3], gd[3];
  
  if (fa) {
    dang = rv3_dih(a, b, c, d, ga, gb, gc, gd) - ang0;
    amp  = (real)( -k1 * sin(dang) - 3 * k3 * sin(3*dang));
    rv3_sinc(fa, ga, amp);
    rv3_sinc(fb, gb, amp);
    rv3_sinc(fc, gc, amp);
    rv3_sinc(fd, gd, amp);
  } else {
    dang = rv3_dih(a, b, c, d, NULL, NULL, NULL, NULL) - ang0;
  }
  return (real)( k1 * (1 - cos(dang)) + k3 * (1 - cos(3 * dang)) );
}

/* 12-10 potential: u = 5(rc/r)^12 - 6(rc/r)^10,
 * the minimum is at r = rc, and u = -1 */
INLINE real pot1210(rv3_t a, rv3_t b, real rc2, real eps, rv3_t fa, rv3_t fb)
{
  real dx[3], dr2, invr2, invr4, invr6, invr10, amp;

  dr2 = rv3_sqr( rv3_diff(dx, a, b) );
  invr2 = rc2 / dr2;
  invr4 = invr2 * invr2;
  invr6 = invr4 * invr2;
  invr10 = invr4 * invr6;
  if (fa) {
    amp = 60 * eps * (invr2 - 1) * invr10 * (1/dr2);
    rv3_sinc(fa, dx,  amp);
    rv3_sinc(fb, dx, -amp);
  }
  return eps * (5 * invr2 - 6) * invr10;
}

/* WCA potential: u = (rc/r)^12 - 2 (rc/r)^6 + 1 if r < rc, or 0 otherwise
 * without the truncation, the minimum is at r = rc, and u = 0 there */
INLINE real potwca(rv3_t a, rv3_t b, real rc2, real eps, rv3_t fa, rv3_t fb)
{
  real dx[3], dr2, invr2, invr6, amp;

  dr2 = rv3_sqr( rv3_diff(dx, a, b) );
  if (dr2 > rc2) return 0;
  invr2 = rc2 / dr2;
  invr6 = invr2 * invr2 * invr2;
  if (fa) {
    amp = 12 * eps * (invr6 - 1) * invr6 * (1/dr2);
    rv3_sinc(fa, dx,  amp);
    rv3_sinc(fb, dx, -amp);
  }
  return eps * (invr6 * (invr6 - 2) + 1);
}

/* repulsive potential: (rc/r)^12 */
INLINE real potr12(rv3_t a, rv3_t b, real rc2, real eps, rv3_t fa, rv3_t fb)
{
  real dx[3], dr2, invr2, invr6, u, amp;

  dr2 = rv3_sqr( rv3_diff(dx, a, b) );
  invr2 = rc2 / dr2;
  invr6 = invr2 * invr2 * invr2;
  u = eps * invr6 * invr6;
  if (fa) {
    amp = 12 * u / dr2;
    rv3_sinc(fa, dx,  amp);
    rv3_sinc(fb, dx, -amp);
  }
  return u;
}


INLINE real cago_force(cago_t *go, rv3_t *x, rv3_t *f)
{
  int i, j, id, n = go->n, ncwca = go->flags & CAGO_NCWCA;
  real ene = 0, kb = go->kb, ka = go->ka, kd1 = go->kd1, kd3 = go->kd3,
       nbe = go->nbe, nbc2 = go->nbc*go->nbc;

  if (f != NULL) {
    for (i = 0; i < n; i++) rv3_zero(f[i]);
  }

  /* bonds */
  for (i = 0; i < n-1; i++)
    ene += potbond(x[i], x[i+1], go->bref[i], kb, f[i], f[i+1]);

  /* angles */
  for (i = 1; i < n-1; i++)
    ene += potang(x[i-1], x[i], x[i+1], go->aref[i-1], ka, f[i-1], f[i], f[i+1]);

  /* dihedrals */
  for (i = 0; i < n - 3; i++)
    ene += potdih13(x[i], x[i+1], x[i+2], x[i+3], go->dref[i],
        kd1, kd3, f[i], f[i+1], f[i+2], f[i+3]);

  /* nonbonded */
  for (i = 0; i < n - 4; i++)
    for (j = i + 4; j < n; j++) {
      id = i*n + j;
      if ( go->iscont[id] ) { /* contact pair */
        ene += pot1210(x[i], x[j], go->r2ref[id], nbe, f[i], f[j]);
      } else { /* noncontact pair */
        ene += ncwca ? potwca(x[i], x[j], nbc2, nbe, f[i], f[j])
          : potr12(x[i], x[j], nbc2, nbe, f[i], f[j]);
      }
    }
  return ene;
}

/* velocity verlet */
INLINE int cago_vv(cago_t *go, real fscal, real dt)
{
  int i, n = go->n;
  real dth = .5f*dt;
  rv3_t *v = go->v, *x = go->x, *f = go->f;

  for (i = 0; i < n; i++) { /* vv part 1 */
    rv3_sinc(v[i], f[i], dth*fscal);
    rv3_sinc(x[i], v[i], dt);
  }

  go->epot = cago_force(go, go->x, go->f); /* calculate force */

  for (i = 0; i < n; i++) { /* vv part 2 */
    rv3_sinc(v[i], f[i], dth*fscal);
  }
  go->ekin = cago_ekin(go, go->v);
  go->t += dt;
  return 0;
}

/* the change of the potential energy */
INLINE real cago_depot(cago_t *go, rv3_t *x, int i, rv3_t xi)
{
  int j, id, n = go->n, ncwca = go->flags & CAGO_NCWCA;
  rv3_t *xn = go->x1; /* we use x1 freely here */
  real ene = 0;
  real ka = go->ka, kb = go->kb, kd1 = go->kd1, kd3 = go->kd3;
  real nbe = go->nbe, nbc2 = go->nbc * go->nbc;
 
  /* copy coordinates */ 
  for (j = 0; j < n; j++) {
    if (j == i) rv3_copy(xn[i], xi);
    else rv3_copy(xn[j], x[j]);
  }

  /* bonds */
  for (j = i - 1; j <= i; j++) {
    if (j < 0 || j >= n - 1) continue;
    ene -= potbond(x[j], x[j+1], go->bref[j], kb, NULL, NULL);
    ene += potbond(xn[j], xn[j+1], go->bref[j], kb, NULL, NULL);
  }

  /* angles */
  for (j = i - 1; j <= i + 1; j++) {
    if (j < 1 || j >= n - 1) continue;
    ene -= potang(x[j-1], x[j], x[j+1], go->aref[j-1], ka,
       NULL, NULL, NULL);
    ene += potang(xn[j-1], xn[j], xn[j+1], go->aref[j-1], ka,
       NULL, NULL, NULL);
  }

  /* dihedrals */
  for (j = i - 3; j <= i; j++) {
    if (j < 0 || j >= n - 3) continue;
    ene -= potdih13(x[j], x[j+1], x[j+2], x[j+3], go->dref[j],
        kd1, kd3, NULL, NULL, NULL, NULL);
    ene += potdih13(xn[j], xn[j+1], xn[j+2], xn[j+3], go->dref[j],
        kd1, kd3, NULL, NULL, NULL, NULL);
  }
  
  /* nonbonded interaction */
  for (j = 0; j < n; j++) {
    if (abs(i - j) < 4) continue;

    /* subtract the old energies */
    id = i*n + j;
    if ( go->iscont[id] ) { /* contact pair */
      ene -= pot1210(x[i], x[j], go->r2ref[id], nbe, NULL, NULL);
    } else { /* noncontact pair */
      ene -= ncwca ? potwca(x[i], x[j], nbc2, nbe, NULL, NULL)
        : potr12(x[i], x[j], nbc2, nbe, NULL, NULL);
    }
    
    /* add the new energies */
    if ( go->iscont[id] ) { /* contact pair */
      ene += pot1210(x[i], x[j], go->r2ref[id], nbe, NULL, NULL);
    } else { /* noncontact pair */
      ene += ncwca ? potwca(x[i], x[j], nbc2, nbe, NULL, NULL)
        : potr12(x[i], x[j], nbc2, nbe, NULL, NULL);
    }
  }
  return ene;
}

/* metropolis algorithm */
INLINE int cago_metro(cago_t *go, real amp, real bet)
{
  int i;
  rv3_t xi;
  real du;

  i = (int) (go->n * rnd0());
  rv3_inc(rv3_rnd(xi, -amp, 2.f*amp), go->x[i]);
  du = cago_depot(go, go->x, i, xi);
  if (du < 0 || rnd0() < exp(-bet * du)) {
    rv3_copy(go->x[i], xi);
    go->epot += du;
    return 1;
  } else {
    return 0;
  }
}

/* write position/velocity file */
INLINE int cago_writepos(cago_t *go, rv3_t *x, rv3_t *v, const char *fn)
{
  FILE *fp;
  int i, n = go->n;

  if (fn == NULL) fn = "cago.pos";
  xfopen(fp, fn, "w", return -1);

  fprintf(fp, "# %d %d\n", go->n, (v != NULL));
  for (i = 0; i < n; i++) {
    fprintf(fp, "%16.12f %16.12f %16.12f ", x[i][0], x[i][1], x[i][2]);
    if (v)
      fprintf(fp, "%16.12f %16.12f %16.12f ", v[i][0], v[i][1], v[i][2]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* count the number of native contacts that are formed in the structure `x'
 * this counting process is independent of the process of defining contacts
 *   although the two processes are very similar, and may be the same
 * here, given a set of defined contacts, we simple observe how many pairs
 *   are close enough to be regarded as contacts
 * a contact is formed if the pair distance is <= gam * native-distance
 * return the number of contacts
 * `*Q' is the ratio of formed contacts / the total number of contacts  */
INLINE int cago_ncontacts(cago_t *go, rv3_t *x, real gam, real *Q, int *mat)
{
  int i, j, id, nct = 0, n = go->n;

  if (gam < 0) gam = 1.2; /* default value */
  if (mat) for (id = 0; id < n * n; id++) mat[id] = 0;

  for (i = 0; i < n - 1; i++)
    for (j = i+1; j < n; j++) {
      if (!go->iscont[ (id = i*n + j) ]) continue;
      if (rv3_dist(x[i], x[j]) < sqrt(go->r2ref[id]) * gam) {
        if (mat) mat[id] = mat[j*n + i] = 1;
        nct++;
      }
    }

  if (Q) *Q = nct / go->ncont;
  return nct;
}

/* read position/velocity file */
INLINE int cago_readpos(cago_t *go, rv3_t *x, rv3_t *v, const char *fn)
{
  char s[1024], *p;
  FILE *fp;
  int i, hasv = 0, next, n = go->n;
  const char *fmt;
  real vtmp[3], *vi;

  if (fn == NULL) fn = "cago.pos";
  xfopen(fp, fn, "r", return -1);

  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "Warning: %s has no information line\n", fn);
    rewind(fp);
  } else {
    if (2 != sscanf(s+1, "%d%d", &i, &hasv) || i != n ) {
      fprintf(stderr, "first line is corrupted:\n%s", s);
      goto ERR;
    }
  }

  if (sizeof(double) == sizeof(real))
    fmt = "%lf%lf%lf%n";
  else
    fmt = "%f%f%f%n";
  for (i = 0; i < n; i++) {
    fgets(s, sizeof s, fp);
    if (strlen(s) < 10) goto ERR;
    if (3 != sscanf(s, fmt, &x[i][0], &x[i][1], &x[i][2], &next)) {
      fprintf(stderr, "cannot read x, i = %d\n", i);
      goto ERR;
    }
    p = s+next;
    if (hasv) {
      vi = (v != NULL) ? v[i] : vtmp;
      if (3 != sscanf(p, fmt, &vi[0], &vi[1], &vi[2], &next)) {
        fprintf(stderr, "cannot read v, i = %d\n", i);
        goto ERR;
      }
    }
  }
  fclose(fp);
  return 0;

ERR:
  fprintf(stderr, "position file [%s] appears to be broken on line %d!\n%s\n", fn, i, s);
  fclose(fp);
  return -1;
}

/* output pdb format */
INLINE int cago_writepdb(cago_t *go, rv3_t *x, const char *fn)
{
  FILE *fp;
  int i, n = go->n;

  xfopen(fp, fn, "w", return -1);
  for (i = 0; i < n; i++)
    fprintf(fp, "ATOM  %5d  CA  %-4sA%4d    %8.3f%8.3f%8.3f  1.00  0.00           C  \n", 
        i+1, pdbaaname(go->aa[i]), i+1, x[i][0], x[i][1], x[i][2]);
  fprintf(fp, "END%77s\n", " ");
  fclose(fp);
  return 0;
}

#endif

