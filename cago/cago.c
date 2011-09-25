#include "util.c"
#include "pdb.c"
#include "rng.c"
#include "dihc.c"
#include "rotfit.c"
#include "md.c"
#include "av.h"
#ifndef CAGO_C__
#define CAGO_C__
#include "cago.h"

/* initialize data for potential energy */
static int cago_initpot(cago_t *go)
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

/* hydro include hydrogen atoms */
static int *pdbcontact(pdbmodel_t *pm, double rc, int hydro)
{
  int ir, jr, i, j, n = pm->n, nres = pm->nres;
  pdbatom_t *at = pm->at;
  real d, dmin;
  int *ds;

  xnew(ds, nres*nres);
  for (ir = 0; ir < nres; ir++) {
    for (jr = ir+1; jr < nres; jr++) {
      /* compute the minimal distance between ir and jr */
      dmin = 1e9;
      for (i = 0; i < n; i++) {
        if (at[i].rid != ir) continue;
        if (!hydro && at[i].atnm[0] == 'H') continue;
        for (j = 0; j < n; j++) {
          if (!hydro && at[j].atnm[0] == 'H') continue;
          if (at[j].rid != jr) continue;
          d = rv3_dist(at[i].x, at[j].x);
          if (d < dmin) dmin = d;
        }
      }
      ds[ir*nres+jr] = ds[jr*nres+ir] = (dmin < rc) ? 1 : 0;
    }
  }

  /* exclude nearby contacts */
  for (ir = 0; ir < nres; ir++)
    for (jr = ir+1; jr < ir+4 && jr < nres; jr++)
      ds[ir*nres + jr] = ds[jr*nres + ir] = 0;
  return ds;
}

/* return cago_t from pdb file fnpdb
 * rcc is the cutoff radius for defining contacts */
cago_t *cago_open(const char *fnpdb, real kb, real ka, real kd1, real kd3,
    real nbe, real nbc, real rcc)
{
  cago_t *go;
  int i, j;
  pdbmodel_t *pm;
  pdbaac_t *c;

  xnew(go, 1);

  if ((pm = pdbm_read(fnpdb, 0)) == NULL)
    return NULL;
  go->iscont = pdbcontact(pm, rcc, 0);
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
  for (i = 0; i < go->n; i++) {
    for (j = 0; j < 3; j++)
      go->xref[i][j] = (real) c->res[i].xca[j];
    go->aa[i] = c->res[i].aa;
  }
  pdbaac_free(c);

  cago_initpot(go);

  { real rmin = 1e9, rmax = 0; int id;
  for (i = 0; i < go->n - 1; i++)
    for (j = i+1; j < go->n; j++)
      if (go->iscont[ (id = i*go->n + j) ]) {
        if (go->r2ref[id] > rmax)
          rmax = go->r2ref[id];
        else if (go->r2ref[id] < rmin)
          rmin = go->r2ref[id];
      }
  rmin = sqrt(rmin); rmax = sqrt(rmax);
  printf("rmsd: %g, %g\n", rmin, rmax); }
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
void cago_rmcom(cago_t *go, rv3_t *x, rv3_t *v)
{
  md_shiftcom3d(x, go->n);
  md_shiftcom3d(v, go->n);
  md_shiftang3d(x, v, go->n);
}

/* initialize a md system */
int cago_initmd(cago_t *go, double rndamp, double T0)
{
  int i, j, n = go->n;
  real s, dx[3];

  xnew(go->x, n);
  xnew(go->v, n);
  xnew(go->f, n);
  xnew(go->x1, n);

  /* initialize position */
  if (rndamp < 0) { /* open chain */
    rndamp *= -1;
    for (i = 0; i < n-1; i++) {
      for (j = 0; j < 3; j++)
        dx[j] = (j == 0) ? 1.f : (rndamp*(2.f*rnd0()/RAND_MAX - 1));
      rv3_normalize(dx);
      rv3_smul(dx, go->bref[i]);
      rv3_add(go->x[i+1], go->x[i], dx);
    }
  } else { /* copy from xref, slightly disturb it */
    for (i = 0; i < n; i++) {
      rv3_copy(go->x[i], go->xref[i]);
      for (j = 0; j < 3; j++)
        go->x[i][j] += rndamp*(2.f*rnd0()/RAND_MAX - 1);
    }
  }
  go->epot = cago_force(go, go->f, go->x);

  /* initialize velocities */
  for (j = 0; j < 3; j++)
    for (i = 0; i < n; i++)
      go->v[i][j] = rnd0() - .5;
  cago_rmcom(go, go->x, go->v); /* remove center of mass motion */
  for (s = 0, i = 0; i < n; i++)
    s += rv3_sqr(go->v[i]);
  s = sqrt( (3*n*T0)/s );
  for (go->ekin = 0, i = 0; i < n; i++) {
    rv3_smul(go->v[i], s);
  }
  go->ekin = cago_ekin(go, go->v);
  go->rmsd = cago_rotfit(go, go->x, NULL);
  go->t = 0;
  //go->mddt = 0.002f;
  //go->thermdt = 0.02f;
  //go->nstcom = 10;
  return 0;
}

real cago_force(cago_t *go, rv3_t *f, rv3_t *x)
{
  int i, j, n = go->n;
  real ene = 0, dr, del, ang, amp, invr2, dr2, dr6;
  real kb = go->kb, ka = go->ka, kd1 = go->kd1, kd3 = go->kd3,
       nbe = go->nbe, nbc2 = go->nbc*go->nbc;
  rv3_t dx, g[3];
  dihcalc_t dc[1];

  for (i = 0; i < n; i++)
    rv3_zero(f[i]);
  /* bonds */
  for (i = 0; i < n-1; i++) {
    rv3_diff(dx, x[i], x[i+1]);
    dr = rv3_norm(dx);
    del = dr - go->bref[i];
    ene += .5f*kb*del*del;
    amp = kb*del/dr;
    rv3_sinc(f[i],   dx, -amp);
    rv3_sinc(f[i+1], dx,  amp);
  }

  /* angles */
  for (i = 1; i < n-1; i++) {
    ang = rv3_ang(x[i-1], x[i], x[i+1], g[0], g[1], g[2]);
    del = ang - go->aref[i-1];
    ene += .5f*ka*del*del;
    amp = -ka*del;
    rv3_sinc(f[i-1], g[0], amp);
    rv3_sinc(f[i],   g[1], amp);
    rv3_sinc(f[i+1], g[2], amp);
  }

  /* dihedrals */
  memset(dc, 0, sizeof(*dc));
  dc->szreal = sizeof(real);
  for (i = 0; i < n - 3; i++) {
    ang = rv3_calcdih(dc, x[i], x[i+1], x[i+2], x[i+3],
                      DIH_FOUR|DIH_GRAD);
    del = ang - go->dref[i];
    if (del > M_PI) del -= 2*M_PI; else if (del < -M_PI) del += 2*M_PI;
    ene += (real)(  kd1 * (1. - cos(del)) );
    amp  = (real)( -kd1 * sin(del) );
    ene += (real)(  kd3 * (1. - cos(3.*del)) );
    amp += (real)( -kd3 * 3 * sin(3*del) );
    rv3_sinc(f[i],   dc->g[0], amp);
    rv3_sinc(f[i+1], dc->g[1], amp);
    rv3_sinc(f[i+2], dc->g[2], amp);
    rv3_sinc(f[i+3], dc->g[3], amp);
  }

  /* nonbonded */
  for (i = 0; i < n - 4; i++) {
    for (j = i + 4; j < n; j++) {
      rv3_diff(dx, x[i], x[j]);
      dr2 = rv3_sqr(dx);
      invr2 = 1/dr2;
      if (go->iscont[i*n + j]) { /* is a contact */
        dr2 = go->r2ref[i*n+j]*invr2;
        dr6 = dr2*dr2*dr2;
        amp = nbe*(60*dr6 - 36)*dr6*invr2;
        ene += nbe*(5*dr6 - 6)*dr6;
      } else {
        dr2 = nbc2/dr2;
        dr6 = dr2*dr2*dr2;
	dr6 *= dr6;
        amp = nbe*12*dr6*invr2;
        ene += nbe*dr6;
      }
      rv3_sinc(f[i], dx,  amp);
      rv3_sinc(f[j], dx, -amp);
    }
  }
  return ene;
}

/* velocity verlet */
int cago_vv(cago_t *go, real fscal, real dt)
{
  int i, n = go->n;
  real dth = .5f*dt;
  rv3_t *v = go->v, *x = go->x, *f = go->f;

  for (i = 0; i < n; i++) { /* vv part 1 */
    rv3_sinc(v[i], f[i], dth*fscal);
    rv3_sinc(x[i], v[i], dt);
  }

  go->epot = cago_force(go, go->f, go->x); /* calculate force */

  for (i = 0; i < n; i++) { /* vv part 2 */
    rv3_sinc(v[i], f[i], dth*fscal);
  }
  go->ekin = cago_ekin(go, go->v);
  go->t += dt;
  return 0;
}

/* write position/velocity file */
int cago_writepos(cago_t *go, rv3_t *x, rv3_t *v, const char *fn)
{
  FILE *fp;
  int i, n = go->n;

  if (fn == NULL) fn = "cago.pos";
  if ((fp = fopen(fn, "w")) == 0) {
    fprintf(stderr, "cannot open file [%s]\n", fn);
    return 1;
  }

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

/* read position/velocity file */
int cago_readpos(cago_t *go, rv3_t *x, rv3_t *v, const char *fn)
{
  char s[1024], *p;
  FILE *fp;
  int i, hasv = 0, next, n = go->n;
  const char *fmt;
  real vtmp[3], *vi;

  if (fn == NULL) fn = "ab.pos";
  if ((fp = fopen(fn, "r")) == 0) {
    fprintf(stderr, "cannot open file [%s]\n", fn);
    return 1;
  }

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
  return 1;
}

/* output pdb format */
int cago_writepdb(cago_t *go, rv3_t *x, const char *fn)
{
  FILE *fp;
  int i, n = go->n;

  if ((fp = fopen(fn, "w")) == 0) {
    fprintf(stderr, "cannot open file [%s]\n", fn);
    return 1;
  }
  for (i = 0; i < n; i++)
    fprintf(fp, "ATOM  %5d  CA  %-4sA%4d    %8.3f%8.3f%8.3f  1.00  0.00           C  \n", 
        i+1, pdbaaname(go->aa[i]), i+1, x[i][0], x[i][1], x[i][2]);
  fprintf(fp, "END%77s\n", " ");
  fclose(fp);
  return 0;
}

/* run a regular md
 * teql steps for equilibration, tmax steps for production
 * tp: the real temperature, tps: thermostat temperature */
int cago_mdrun(cago_t *go, real mddt, real thermdt, int nstcom, 
    real tps, real tp, av_t *avep, av_t *avrmsd,
    int teql, int tmax, int trep)
{
  int t;
  real fs = tps/tp;

  tmax = (tmax < 0) ? -1 : (tmax + teql);
  av_clear(avep);
  av_clear(avrmsd);
  for (t = 1; tmax < 0 || t <= tmax; t++) {
    cago_vv(go, fs, mddt);
    if (t % nstcom == 0) cago_rmcom(go, go->x, go->v);
    cago_vrescale(go, (real) tps, thermdt);
    go->rmsd = cago_rotfit(go, go->x, NULL);
    if (t > teql) {
      av_add(avep, go->epot);
      av_add(avrmsd, go->rmsd);
    }
    if (trep > 0 && t % trep == 0) {
      printf("%9d: tp %.4f, tps %.4f, rmsd %7.4f, K %.2f, U %.2f\n",
          t, tp, tps, go->rmsd, go->ekin, go->epot);
    }
  }
  return 0;
}

/* guess a proper temperature for a given rmsd,
 * return 0 if successful.
 *
 * temperature is updated according to rmsd 
 * several stages of updating are used, each with a fixed tpdt
 * after a stage, the updating magnitude amp is multiplied by ampf
 * iterations finish when the temperature difference is less than 
 * a given tolerance 'tptol'
 * a pass is defined every time the rmsd crosses 'rmsd'
 * in every stage, npass passes are required to determine convergence
 * */
int cago_cvgmdrun(cago_t *go, real mddt, real thermdt, int nstcom,
    real rmsd, int npass, 
    real amp, real ampf, real tptol, av_t *avtp, av_t *avep,
    real tp, real tpmin, real tpmax, int tmax, int trep)
{
  int i, t, stg, sgp, sgn, ipass;
  real tpp = 0, tp1, tpav, tmp;

  go->rmsd = cago_rotfit(go, go->x, NULL);
  sgp = (go->rmsd > rmsd) ? 1 : -1;
  for (stg = 0; ; stg++, amp *= ampf) { /* stages with different dpdt */
    if (avtp) av_clear(avtp);
    if (avep) av_clear(avep);
    for (ipass = 0, t = 1; (tmax < 0 || t <= tmax) && ipass < npass; t++) {
      cago_vv(go, 1, mddt);
      if (t % nstcom == 0) cago_rmcom(go, go->x, go->v);
      cago_vrescale(go, (real) tp, thermdt);
      go->rmsd = cago_rotfit(go, go->x, NULL);
      sgn = (go->rmsd > rmsd) ? 1 : -1;
      if (sgn * sgp < 0) {
        ipass++;
        sgp = sgn;
      }
      /* update the temperature */
      tp1 = tp - sgn*mddt*amp;
      if (tp1 < tpmin) tp1 = tpmin;
      else if (tp1 > tpmax) tp1 = tpmax;
      for (tmp = tp1/tp, i = 0; i < go->n; i++) /* scale v */
        rv3_smul(go->v[i], tmp);
      tp = tp1;
      if (avtp) av_add(avtp, tp);
      if (avep) av_add(avep, go->epot);
      if (trep >= 0 && t % trep == 0) {
        printf("%d|%9d: %.2f - %.2f, tp %.4f, K %.2f, U %.2f, pass: %d/%d\n",
            stg, t, go->rmsd, rmsd, tp,
            go->ekin, go->epot, ipass, npass);
      }
    }
    /* end of a stage */
    if (ipass < npass) { /* not enough passes over rmsd */
      const char fnfail[] = "fail.pos";
      cago_rotfit(go, go->x, go->x1);
      cago_writepos(go, go->x1, NULL, fnfail);
      fprintf(stderr, "%d: failed to converge, rmsd: %g - %g, %d passes, %s\n", 
          stg, rmsd, go->rmsd, ipass, fnfail);
      return 1;
    }
    tpav = av_getave(avtp);
    printf("%d: amp %g, tp %g, tpav %g/%g, epotav %g, pass %d/%d\n", 
        stg, amp, tp, tpav, tpp, av_getave(avep), ipass, npass);
    tmp = .5*(tpav + tpp);
    if (stg > 0 && fabs(tpav - tpp) < tptol*tmp) break;
    tpp = tpav;
  }
  return 0;
}



#endif

