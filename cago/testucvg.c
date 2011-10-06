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

real epot_target = 0;
int npass = 400;

static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "C-alpha GO model potential energy convergent";
  argopt_regarg(ao, NULL, &fnpdb, "pdbfile");
  argopt_regopt_help(ao, "-h");
  argopt_regopt(ao, "-E", "%r", &epot_target, "target energy");
  argopt_regopt(ao, "-p","%d", &epot_target, "number of passes");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}

int main(int argc, char **argv)
{
  cago_t *go;
  int ret;
  int nstcom = 10, tmax = 100000000, trep = 10000;
  real tptol = 0.01f, amp = 0.01f, ampf = sqrt(0.1);
  real mddt = 0.002f, thermdt = 0.02f;
  real tptry = 1.0, tpmin = 0.001, tpmax = 50.0;
  av_t avtp[1], avep[1], avrmsd[1];
  real tpav, tpdv, epav, epdv, rdav, rddv;

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, 0.0);
  printf("%s n %d, epot = %g, %g, rmsd = %g\n", fnpdb, go->n, go->epot, go->epotref, go->rmsd);
  if (epot_target < go->epotref) {
    fprintf(stderr, "target energy %g must be greater than %g\n", epot_target, go->epotref);
    return 0;
  }

  ret = cago_ucvgmdrun(go, mddt, thermdt, nstcom,
      epot_target, npass, amp, ampf, tptol, avtp, avep, avrmsd, 
      tptry, tpmin, tpmax, tmax, trep);
  tpav = av_getave(avtp);
  tpdv = av_getdev(avtp);
  epav = av_getave(avep);
  epdv = av_getdev(avep);
  rdav = av_getave(avrmsd);
  rddv = av_getdev(avrmsd);
  printf("ucvgmd %s, epot %7.4f, tp = %.4f(%.4f), epot = %.4f(%.4f), rmsd = %.2f(%.2f)\n",
      (ret ? "   failed" : "succeeded"), epot_target, tpav, tpdv, epav, epdv, rdav, rddv);
  cago_rotfit(go, go->x, go->f);
  cago_writepos(go, go->f, NULL, "c.pos");
  cago_writepos(go, go->xref, NULL, "ref.pos");  
  cago_close(go);
  return 0;
}
