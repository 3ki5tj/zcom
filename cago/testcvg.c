#include "cago.c"

const char *fnpdb = "pdb/1MBN.pdb";
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = .5f;
real nbe = 1.f;
real nbc = 4.f; /* repulsion distance */
real rcutoff = 6.f;

real tp = 1.0f;
real mddt = 2e-3f;
real thermdt = 2e-2f;

int tfreq =  2000;
int tmax = 200000;

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
  int ret, npass = 10;
  double rmsd_target = 15.0, tptol = 0.1, amp = 0.01, ampf = sqrt(0.1);
  av_t avtp[1], avep[1];

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcutoff)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, 0.0);
  printf("ene = %g, %g, rmsd = %g\n", go->epot, go->ekin, go->rmsd);
  ret = cago_cvgmdrun(go, rmsd_target, npass, 
      amp, ampf, tptol, avtp, avep, 
      1.0, 0.01, 100.0, 1000000, 1000);
  printf("cvgmd %s, epot = %g, tp = %g\n", 
      (ret ? "failed" : "succeeded"), 
      av_getave(avep), av_getave(avtp));
  cago_close(go);
  return 0;
}
