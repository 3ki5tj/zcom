#include "cago.c"

const char *fnpdb = "pdb/1VII.pdb";
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = .5f;
real nbe = 1.f;
real nbc = 4.f; /* repulsion distance */
real rcc = 6.f;

const char *prog = "cagocvg";

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
