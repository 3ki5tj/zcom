#include "cago.c"
#include "argopt.c"

const char *fnpdb = "pdb/1KIK.pdb";
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = .5f;
real nbe = 1.f;
real nbc = 4.f; /* repulsion distance */
real rcc = 4.5f;

real tps = 1.12f, tp = 1.12f;

static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  argopt_regarg(ao, NULL, &fnpdb, "pdbfile");
  argopt_reghelp(ao, "-h");
  argopt_regopt(ao, "-T", "%r", &tp, "Temperature");
  argopt_parse(ao, argc, argv);
  if (argopt_set(ao, tp)) tps = tp;
  argopt_close(ao);
}

/* run a regular md
 * teql steps for equilibration, tmax steps for production
 * tp: the real temperature, tps: thermostat temperature */
INLINE int cago_mdrun(cago_t *go, real mddt, real thermdt,
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
    cago_rmcom(go, go->x, go->v);
    cago_vrescale(go, (real) tps, thermdt);
    go->rmsd = cago_rmsd(go, go->x, NULL);
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

int main(int argc, char **argv)
{
  cago_t *go;
  int teql = 100000, tmax = 500000, trep = 10000, ncont;
  real mddt = 0.002f, thermdt = 0.02f;
  real epav, epdv, rdav, rddv;
  av_t avep[1], avrmsd[1];

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, -0.1, 0.0);
  ncont = cago_ncontacts(go, go->x, 1.2, NULL, NULL);
  printf("%s n %d, tp %.3f, tps %.3f, epot %g, %g (ref), rmsd %g, cont %d/%d\n", 
      fnpdb, go->n, tp, tps, go->epot, go->epotref,
      go->rmsd, ncont, go->ncont);
  
  cago_mdrun(go, mddt, thermdt, tps, tp, avep, avrmsd,
     teql, tmax, trep);
  epav = av_getave(avep);
  epdv = av_getdev(avep);
  rdav = av_getave(avrmsd);
  rddv = av_getdev(avrmsd);
  printf("tp %.3f, tps %.3f, epot %.2f(%.2f), rmsd %.4f(%.4f)\n", 
      tp, tps, epav, epdv, rdav, rddv);
  cago_close(go);
  return 0;
}
