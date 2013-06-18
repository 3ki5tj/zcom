#include "av.h"
#include "argopt.h"
#include "cago.h"


const char *fnpdb = "pdb/1VII.pdb";
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = .5f;
real nbe = 1.f;
real nbc = 4.f; /* repulsion distance */
real rcc = 6.f;

int dormsd = 0;
real epottarget = 0;
real rmsdtarget = 5;  /* angstroms */
int npass = 100; /* number of times to cross over the boundary
                    and back to determine convergence */


static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "C-alpha Go model";
  argopt_add(ao, NULL, NULL, &fnpdb, "pdbfile");
  argopt_addhelp(ao, "-h");
  argopt_add(ao, "-r", "%b", &dormsd, "do RMSD instead of potential energy");
  argopt_add(ao, "-E", "%r", &epottarget, "target energy");
  argopt_add(ao, "-R", "%r", &rmsdtarget, "target RMSD");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}



/* guess a proper temperature for a target potential energy or
 * a target RMSD
 * return 0 if successful
 *
 * temperature is updated according to epot
 * several stages of updating are used, each with a fixed tpdt
 * after a stage, the updating magnitude amp is multiplied by ampf
 * iterations finish when the temperature difference is less than
 * a given tolerance 'tptol'
 * a pass is defined every time the potential energy crosses 'epot'
 * in every stage, npass passes are required to determine convergence
 * */
INLINE int cago_cvgmdrun(cago_t *go, real mddt, real thermdt, int nstcom,
    int dormsd, real epot, real rmsd, int npass,
    real amp, real ampf, real tptol, av_t *avtp, av_t *avep, av_t *avrmsd,
    real tp, real tpmin, real tpmax, int tmax, int trep)
{
  int i, t, stg, sgp, sgn, ipass;
  real tpp = 0, tp1, tpav, epav, rdav, tmp;

  go->rmsd = cago_rmsd(go, go->x, NULL);
  if (dormsd) sgp = (go->rmsd > rmsd) ? 1 : - 1;
  else sgp = (go->epot > epot) ? 1 : -1;
  for (stg = 0; ; stg++, amp *= ampf) { /* stages with different dpdt */
    if (avtp) av_clear(avtp);
    if (avep) av_clear(avep);
    if (avrmsd) av_clear(avrmsd);
    for (ipass = 0, t = 1; (tmax < 0 || t <= tmax) && ipass < npass; t++) {
      cago_vv(go, 1, mddt);
      if (t % nstcom == 0) cago_rmcom(go, go->x, go->v);
      cago_vrescalex(go, (real) tp, thermdt, &go->ekin, &go->tkin);
      go->rmsd = cago_rmsd(go, go->x, NULL);
      if (dormsd) sgn = (go->rmsd > rmsd) ? 1 : -1;
      else sgn = (go->epot > epot) ? 1 : -1;
      if (sgn * sgp < 0) {
        ipass++;
        sgp = sgn;
      }

      /* update the temperature by negative feedback */
      tp1 = tp - sgn * mddt * amp;

      if (tp1 < tpmin) tp1 = tpmin;
      else if (tp1 > tpmax) tp1 = tpmax;
      for (tmp = tp1/tp, i = 0; i < go->n; i++) /* scale v */
        rv3_smul(go->v[i], tmp);
      tp = tp1;
      if (avtp) av_add(avtp, tp);
      if (avep) av_add(avep, go->epot);
      if (avrmsd) av_add(avrmsd, go->rmsd);
      if (trep >= 0 && t % trep == 0) {
        printf("%s %d|%9d: U %.2f - %.2f, R %.1f - %.1f, "
            "tp %.4f, K %.2f, pass: %d/%d\n",
            dormsd ? "R" : "U", stg, t, go->epot, epot,
            go->rmsd, rmsd, tp, go->ekin, ipass, npass);
      }
    }
    /* end of a stage */
    if (ipass < npass) { /* not enough passes over rmsd */
      const char fnfail[] = "fail.pos";
      go->rmsd = cago_rmsd(go, go->x, go->x1);
      cago_writepos(go, go->x1, NULL, fnfail);
      fprintf(stderr, "%d: failed to converge, epot: %g - %g, %g - %g, "
          "%d passes, %s\n",
          stg, epot, go->epot, rmsd, go->rmsd, ipass, fnfail);
      return 1;
    }
    tpav = av_getave(avtp);
    epav = av_getave(avep);
    rdav = av_getave(avrmsd);
    printf("%d: amp %g, tp %g, tpav %g/%g, epotav %g, rmsdav %g, pass %d/%d\n",
        stg, amp, tp, tpav, tpp, epav, rdav, ipass, npass);
    tmp = .5*(tpav + tpp);
    if (stg > 0 && fabs(tpav - tpp) < tptol*tmp) break;
    tpp = tpav;
  }
  return 0;
}



int main(int argc, char **argv)
{
  cago_t *go;
  int ret;
  int nstcom = 10, tmax = 100000000, trep = 10000;
  real tptol = 0.01f, amp = 0.001f, ampf = sqrt(0.1);
  real mddt = 0.002f, thermdt = 0.02f;
  real tp0 = 1.0, tpmin = 0.1, tpmax = 10.0;
  av_t avtp[1], avep[1], avrmsd[1];
  real tpav, tpdv, epav, epdv, rdav, rddv;

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc,
                      PDB_CONTACT_HEAVY, 4, CAGO_VERBOSE)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, 0.0);
  printf("%s n %d, epot = %g, %g, rmsd = %g\n", fnpdb, go->n, go->epot, go->epotref, go->rmsd);
  if (epottarget < go->epotref) {
    fprintf(stderr, "target energy %g must be greater than %g\n", epottarget, go->epotref);
    return 0;
  }

  ret = cago_cvgmdrun(go, mddt, thermdt, nstcom,
      dormsd, epottarget, rmsdtarget,
      npass, amp, ampf, tptol, avtp, avep, avrmsd,
      tp0, tpmin, tpmax, tmax, trep);
  tpav = av_getave(avtp);
  tpdv = av_getdev(avtp);
  epav = av_getave(avep);
  epdv = av_getdev(avep);
  rdav = av_getave(avrmsd);
  rddv = av_getdev(avrmsd);
  printf("%s cvgmd %s epot_t %.4f, rmsd_t %.4f, tp = %.4f(%.4f), "
         "epot %.4f(%.4f), rmsd %.2f(%.2f)\n",
         (dormsd ? "rmsd" : "epot"), (ret ? "failed" : "succeeded"),
         epottarget, rmsdtarget, tpav, tpdv, epav, epdv, rdav, rddv);
  cago_close(go);
  return 0;
}

