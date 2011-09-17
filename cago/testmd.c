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
  int t;

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcutoff)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, 0.1, 0.0);
  printf("ene = %g, %g, rmsd = %g\n", go->epot, go->ekin, go->rmsd);
  cago_writepos(go, go->x, go->v, "a.pos");
  for (t = 1; t <= tmax; t++) {
    tp = (2 - 1.9f*t/tmax);
    cago_vv(go, 1.f, mddt);
    cago_vrescale(go, tp, thermdt);
    cago_rmcom(go, go->x, go->v);
    if (t % tfreq == 0) {
      cago_rotfit(go, go->x, NULL);
      printf("t %d, tp = %g, ene = %g+%g = %g, %g\n", 
    	  t, tp, go->epot, go->ekin, go->epot + go->ekin, go->rmsd);
    }
  }

  cago_rotfit(go, go->x, go->f);
  cago_writepos(go, go->x, NULL, "b.pos");
  cago_writepos(go, go->f, NULL, "c.pos");
  cago_writepos(go, go->xref, NULL, "0.pos");
  cago_writepdb(go, go->f, "final.pdb");
  cago_writepdb(go, go->xref, "ref.pdb");
  printf("ene = %g, %g, rmsd = %g\n", go->epot, go->ekin, go->rmsd);

  cago_close(go);

  return 0;
}
