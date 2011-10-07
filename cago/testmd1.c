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

real tps = 0.3f, tp = 0.3f;

static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  argopt_regarg(ao, NULL, &fnpdb, "pdbfile");
  argopt_regopt_help(ao, "-h");
  argopt_regopt(ao, "-T", "%r", &tp, "Temperature");
  argopt_parse(ao, argc, argv);
  if (argopt_set(ao, tp)) tps = tp;
  argopt_close(ao);
}

int main(int argc, char **argv)
{
  cago_t *go;
  int nstcom = 10, teql = 100000, tmax = 500000, trep = 10000;
  real mddt = 0.002f, thermdt = 0.02f;
  real epav, epdv, rdav, rddv;
  av_t avep[1], avrmsd[1];

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, 0.0);
  printf("%s n %d, tp %.3f, tps %.3f, epot = %g, %g (ref), rmsd = %g\n", 
      fnpdb, go->n, tp, tps, go->epot, go->epotref, go->rmsd);
  
  cago_mdrun(go, mddt, thermdt, nstcom, tps, tp, avep, avrmsd,
     teql, tmax, trep);
  epav = av_getave(avep);
  epdv = av_getdev(avep);
  rdav = av_getave(avrmsd);
  rddv = av_getdev(avrmsd);
  printf("tp %.3f, tps %.3f, epot %.2f(%.2f), rmsd %.4f(%.4f)\n", 
      tp, tps, epav, epdv, rdav, rddv);
  cago_rotfit(go, go->x, go->f);
  cago_writepos(go, go->f, NULL, "c.pos");
  cago_writepos(go, go->xref, NULL, "ref.pos");

  cago_close(go);

  return 0;
}
