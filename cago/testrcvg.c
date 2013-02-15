#include "cago.c"
#include "argopt.c"

const char *fnpdb = "pdb/1VII.pdb";
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = .5f;
real nbe = 1.f;
real nbc = 4.f; /* repulsion distance */
real rcc = 6.f;

static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "C-alpha Go model, rmsd convergent";
  argopt_regarg(ao, NULL, &fnpdb, "pdbfile");
  argopt_reghelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}

/* guess a proper temperature for a target rmsd, return 0 if successful
 * should work for a small rmsd
 * for a rmsd in the transition region, it can be stuck in a local minimal,
 * in which rmsd is greater than the target value, but tp reaches tpmin
 *
 * temperature is updated according to rmsd 
 * several stages of updating are used, each with a fixed tpdt
 * after a stage, the updating magnitude amp is multiplied by ampf
 * iterations finish when the temperature difference is less than 
 * a given tolerance 'tptol'
 * a pass is defined every time the rmsd crosses 'rmsd'
 * in every stage, npass passes are required to determine convergence
 * */
INLINE int cago_rcvgmdrun(cago_t *go, real mddt, real thermdt, int nstcom,
    real rmsd, int npass, 
    real amp, real ampf, real tptol, av_t *avtp, av_t *avep, av_t *avrmsd,
    real tp, real tpmin, real tpmax, int tmax, int trep)
{
  int i, t, stg, sgp, sgn, ipass;
  real tpp = 0, tp1, tpav, epav, rdav, tmp;

  go->rmsd = cago_rmsd(go, go->x, NULL);
  sgp = (go->rmsd > rmsd) ? 1 : -1;
  for (stg = 0; ; stg++, amp *= ampf) { /* stages with different dpdt */
    if (avtp) av_clear(avtp);
    if (avep) av_clear(avep);
    if (avrmsd) av_clear(avrmsd);
    for (ipass = 0, t = 1; (tmax < 0 || t <= tmax) && ipass < npass; t++) {
      cago_vv(go, 1, mddt);
      if (t % nstcom == 0) cago_rmcom(go, go->x, go->v);
      cago_vrescale(go, (real) tp, thermdt);
      go->rmsd = cago_rmsd(go, go->x, NULL);
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
      if (avrmsd) av_add(avrmsd, go->rmsd);
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
    epav = av_getave(avep);
    rdav = av_getave(avrmsd);
    printf("%d: amp %g, tp %g, tpav %g/%g, epotav %g, rmsd %g, pass %d/%d\n", 
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
  int ret, npass = 400;
  int nstcom = 10, tmax = 20000000, trep = 10000;
  real rmsd_target = 0.5f, tptol = 0.01f, amp = 0.01f, ampf = sqrt(0.1);
  real mddt = 0.002f, thermdt = 0.02f;
  real tptry = 0.1, tpmin = 0.01, tpmax = 50.0;
  av_t avtp[1], avep[1], avrmsd[1];
  real tpav, tpdv, epav, epdv, rdav, rddv;

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, 0.0);
  printf("%s n %d, epot = %g, %g, rmsd = %g\n", fnpdb, go->n, go->epot, go->epotref, go->rmsd);

  ret = cago_rcvgmdrun(go, mddt, thermdt, nstcom,
      rmsd_target, npass, amp, ampf, tptol, avtp, avep, avrmsd, 
      tptry, tpmin, tpmax, tmax, trep);
  tpav = av_getave(avtp);
  tpdv = av_getdev(avtp);
  epav = av_getave(avep);
  epdv = av_getdev(avep);
  rdav = av_getave(avrmsd);
  rddv = av_getdev(avrmsd);
  printf("rcvgmd %s, rmsd %7.4f, tp = %.4f(%.4f), epot = %.4f(%.4f), rmsd = %.2f(%.2f)\n",
      (ret ? "   failed" : "succeeded"), rmsd_target, tpav, tpdv, epav, epdv, rdav, rddv);
  cago_rotfit(go, go->x, go->f);
  cago_writepos(go, go->f, NULL, "c.pos");
  cago_writepos(go, go->xref, NULL, "ref.pos");  
  cago_close(go);
  return 0;
}
