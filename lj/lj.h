#include "def.h"
#include "util.h"
#include "hist.h"
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
  double t;

  int usesw; /* switched potential */
  real rs, a4, a5, a6, a7; /* parameters */
  ljpair_t *pr;
  int npr;
  real lap, f2, *gdg, *xdg;

  int usesq; /* square well potential */
  int esqinf;
  real ra, ra2, rb, rb2; /* -1 for (ra, rb) */
  real rmin; /* minimal pair distance */

  hist_t *rdf; /* histogram for radial distribution function */
  int rdfnfr; /* number of frames in rdf */

  unsigned isclone; /* is a clone copy, don't free pointers */
} lj_t;

/* copy flags */
#define LJ_CPX   0x0001
#define LJ_CPV   0x0002
#define LJ_CPF   0x0004
#define LJ_CPRDF 0x0010
#define LJ_CPPR  0x0020
#define LJ_CPGDG 0x0040
#define LJ_CPXDG 0x0080
#define LJ_CPXVF (LJ_CPX|LJ_CPV|LJ_CPF)

lj_t *lj_open(int n, int d, real rho, real rcdef);
void lj_close(lj_t *lj);
INLINE lj_t *lj_copy(lj_t *dest, const lj_t *src, unsigned flags);
INLINE lj_t *lj_clone(const lj_t *src, unsigned flags);

INLINE int lj_writepos(lj_t *lj, const real *x, const real *v, const char *fn);
#define LJ_LOADBOX 0x10
INLINE int lj_readpos(lj_t *lj, real *x, real *v, const char *fn, unsigned flags);

/* open rdf */
INLINE hist_t *lj_rdfopen(lj_t *lj, double dr, double rmax);
/* add pairs to the RDF data */
INLINE int lj_rdfadd(lj_t *lj);
/* save rdf, flags can have HIST_NOZEROES */
INLINE int lj_rdfsave(lj_t *lj, const char *fn, unsigned flags);
/* load rdf, flags can have HIST_ADDITION and/or HIST_VERBOSE */
INLINE int lj_rdfload(lj_t *lj, const char *fn, unsigned flags);

/* compute the equation of states */
INLINE double lj_eos3d(double rho, double T, double *P, double *Ar, double *Gr);

/* initialize the square well potential */
INLINE void lj_initsq(lj_t *lj, real ra, real rb);

/* initialize the switched potential */
INLINE void lj_initsw(lj_t *lj, real rs);

real lj_energy(lj_t *lj);

INLINE real lj_energyx(lj_t *lj, real *x, real *vir, int *iep, real *rmin,
    real *ep0, real *eps, real *lap);

real lj_force(lj_t *lj);

/* Metropolis algorithm */
INLINE int lj_metro(lj_t *lj, real amp, real bet);

/* returns the energy change of a local perturbation */
INLINE real lj_dupertl(lj_t *lj, real amp);

/* return the energy change by an all-atom perturbation */
INLINE real lj_dupertg(lj_t *lj, real amp);

/* compute the configuration temperature */
INLINE real lj_bconfsw3d(lj_t *lj, real *udb);

/* compute pressure */
INLINE real lj_calcp(lj_t *lj, real tp)
  { return (lj->dof * tp + lj->vir) / (lj->d * lj->vol) + lj->p_tail; }

/* compute pressure, ideal gas part from the kinetic energy  */
INLINE real lj_calcpk(lj_t *lj)
  { return (2.f * lj->ekin + lj->vir) / (lj->d * lj->vol) + lj->p_tail; }

/* set density, compute tail corrections, etc */
INLINE void lj_setrho(lj_t *lj, real rho);

#define lj_vv(lj, dt) lj_vvx(lj, dt, 1.f)
INLINE void lj_vvx(lj_t *lj, real dt, real fscal);

#define lj_shiftcom(lj, v)    md_shiftcom(v, lj->n, lj->d)
#define lj_shiftang(lj, x, v) md_shiftang(x, v, lj->n, lj->d)

INLINE void lj_vscale(lj_t *lj, real tp, real ekt) 
 { md_vscale(lj->v, lj->n * lj->d, lj->dof, tp, ekt, &lj->ekin, &lj->tkin); }

INLINE void lj_vrescale(lj_t *lj, real tp, real thermdt) 
 { md_vrescale(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin); }

INLINE void lj_vrescalex(lj_t *lj, real tp, real thermdt)
 { md_vrescalex(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin); }

INLINE int lj_mcvrescale(lj_t *lj, real tp, real thermdt)
 { return md_mcvrescale(lj->v, lj->n * lj->d, lj->dof, tp, thermdt, &lj->ekin, &lj->tkin); }

/* Nose-Hoover thermostat/barostat
 * set cutoff to half of the box */
INLINE void lj_hoovertp(lj_t *lj, real dt, real tp, real pext,
    real *zeta, real *eta, real Q, real W, int ensx)
{
  md_hoovertp(lj->v, lj->n, lj->d, lj->dof, dt, tp, pext, zeta, eta,
      Q, W, lj->vol, lj->vir, lj->p_tail, ensx, &lj->ekin, &lj->tkin);
}

/* Nose-Hoover chain thermostat/barostat
 * set cutoff to half of the box */
INLINE void lj_nhchaintp(lj_t *lj, real dt, real tp, real pext,
    real *zeta, real *eta, const real *Q, int M, real W, int ensx)
{
  md_nhchaintp(lj->v, lj->n, lj->d, lj->dof, dt, tp, pext, zeta, eta,
      Q, M, W, lj->vol, lj->vir, lj->p_tail, ensx, &lj->ekin, &lj->tkin);
}

/* velocity verlet with the scaling step in the Nose-Hoover barostat */
INLINE void lj_vv_hoovertp(lj_t *lj, real dt, real eta);

/* Berendsen barostat: as a backup for a constant pressure simulation */
INLINE void lj_pberendsen(lj_t *lj, real barodt, real tp, real pext);

/* Langevin barostat, with kinetic-energy scaling */
INLINE void lj_langtp(lj_t *lj, real dt, real tp, real pext,
    real zeta, real *eta, real W, int ensx)
{
  md_langtp(lj->v, lj->n, lj->d, dt, tp, pext, zeta, eta,
      W, lj->vol, lj->vir, lj->p_tail, ensx, &lj->ekin, &lj->tkin);
}

/* position Langevin barostat, with kinetic-energy scaling */
INLINE void lj_langtp0(lj_t *lj, real barodt, real tp, real pext, int ensx)
{
  md_langtp0(lj->v, lj->n, lj->d, barodt, tp, pext, &lj->vol,
     lj->vir, lj->p_tail, ensx, &lj->ekin, &lj->tkin);
  lj_setrho(lj, lj->n/lj->vol);
  lj_force(lj);
}

/* old interface */
#define lj_lgvvolmove(lj, barodt, tp, p) lj_langp0(lj, barodt, tp, p, 0)

/* Langevin barostat, with coordinates only, barodt ~ 1e-5 for n = 108 */
INLINE void lj_langp0(lj_t *lj, real barodt, real tp, real pext, int ensx)
{
  md_langp0(lj->dof, lj->d, barodt, tp, pext, &lj->vol, lj->vir, lj->p_tail, ensx);
  lj_setrho(lj, lj->n/lj->vol);
  lj_force(lj);
}

/* In Monte Carlo barostats, we compute the energy directly */
#define LJ_FIXEDRC 0x4000

#define lj_mcprescale(lj, lnvamp, tp, pext, vmin, vmax, ensx) \
  lj_mctp(lj, lnvamp, tp, pext, vmin, vmax, ensx, 0)

/* Monte Carlo barostat, with kinetic-energy scaling */
INLINE int lj_mctp(lj_t *lj, real lnvamp, real tp, real pext,
    real vmin, real vmax, int ensx, unsigned flags);

/* old interface */
#define lj_volmove(lj, lnlamp, tp, p) \
  lj_mcp(lj, lnlamp*lj->d, tp, p, 0, 1e300, 0, LJ_FIXEDRC)

/* Monte Carlo barostat, coordinate only */
INLINE int lj_mcp(lj_t *lj, real vamp, real tp, real pext,
    real vmin, real vmax, int ensx, unsigned flags);

#endif

