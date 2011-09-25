#include "cago.c"

const char *fnpdb = "pdb/1VII.pdb";
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = .5f;
real nbe = 1.f;
real nbc = 4.f; /* repulsion distance */
real rcc = 6.f;

const char *prog = "go";

static void help(void)
{
  printf("%s your.pdb\n", prog);
  exit(1);
}

static void doargs(int argc, const char **argv)
{
  int i;
  
  prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      fnpdb = argv[i];
      continue;
    }
    if (argv[i][1] == 'h')
      help();
  }
}

int main(int argc, const char **argv)
{
  cago_t *go;
  int nstcom = 10, teql = 100000, tmax = 5000000, trep = 10000;
  real mddt = 0.002f, thermdt = 0.02f, tps = 1.1f, tp = 1.1f;
  real epav, epdv, rdav, rddv;
  av_t avep[1], avrmsd[1];

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, 0.0);
  printf("tp %.2f, tps %.2f, epot = %g, %g, rmsd = %g\n", 
      tp, tps, go->epot, go->epotref, go->rmsd);
  
  cago_mdrun(go, mddt, thermdt, nstcom, tps, tp, avep, avrmsd,
     teql, tmax, trep);
  epav = av_getave(avep);
  epdv = av_getdev(avep);
  rdav = av_getave(avrmsd);
  rddv = av_getdev(avrmsd);
  printf("tp %2.f, tps %.2f, epot %.2f(%.2f), rmsd %.4f(%.4f)\n", 
      tp, tps, epav, epdv, rdav, rddv);
  cago_rotfit(go, go->x, go->f);
  cago_writepos(go, go->f, NULL, "c.pos");
  cago_writepos(go, go->xref, NULL, "ref.pos");

  cago_close(go);

  return 0;
}
