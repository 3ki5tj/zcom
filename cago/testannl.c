#include "argopt.h"
#include "cago.h"

const char *fnpdb = "pdb/1VII.pdb";
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = .5f;
real nbe = 1.f;
real nbc = 4.f; /* repulsion distance */
real rc = 5.f; /* cutoff distance of defining contacts */

real tp = 1.0f;
real mddt = 2e-3f;
real thermdt = 0.1f;

int tfreq = 2000;
int tmax = 20000;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, NULL, NULL, &fnpdb, "pdbfile");
  argopt_addhelp(ao, "-h");
  argopt_add(ao, "-n", "%d", &tmax, "number of simulation steps");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}



/* test if energy and force matches with each other */
static void eftest(cago_t *go, real del)
{
  real ep0, ep1, dep, f2;
  int i;

  ep0 = cago_force(go, go->x, go->f);
  for (f2 = 0, i = 0; i < go->n; i++) f2 += rv3_sqr(go->f[i]);
  /* displace the molecule */
  for (i = 0; i < go->n; i++) {
    rv3_sinc(go->x[i], go->f[i], del/f2);
  }
  ep1 = cago_force(go, go->x, go->f);
  dep = ep0 - ep1;
  printf("ep %g, %g, %g, rat %g\n", ep0, ep1, dep, dep/del);
}

int main(int argc, char **argv)
{
  cago_t *go;
  int t;

  doargs(argc, argv);
  if ((go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rc,
                      PDB_CONTACT_HEAVY, 4, CAGO_VERBOSE)) == NULL) {
    fprintf(stderr, "cannot initialize from %s\n", fnpdb);
    return 1;
  }
  cago_initmd(go, -0.1, 0.0);
  printf("ene = %g, %g, rmsd = %g\n", go->epot, go->ekin, go->rmsd);
  cago_writepos(go, go->x, go->v, "a.pos");
  for (t = 1; t <= tmax; t++) {
    tp = (2 - 1.9f*t/tmax); /* annealing */
    cago_vv(go, 1.f, mddt);
    cago_rmcom(go, go->x, go->v);
    cago_vrescale(go, tp, thermdt, &go->ekin, &go->tkin);
    if (t % tfreq == 0) {
      eftest(go, 0.1);
      go->rmsd = cago_rmsd(go, go->x, NULL);
      printf("t %d, tp = %g, ene = %g+%g = %g, %g\n",
    	  t, tp, go->epot, go->ekin, go->epot + go->ekin, go->rmsd);
    }
  }

  go->rmsd = cago_rmsd(go, go->x, go->f);
  cago_writepos(go, go->x, NULL, "b.pos");
  cago_writepos(go, go->f, NULL, "c.pos");
  cago_writepos(go, go->xref, NULL, "0.pos");
  cago_writepdb(go, go->f, "final.pdb");
  cago_writepdb(go, go->xref, "ref.pdb");
  printf("ene = %g, %g, rmsd = %g\n", go->epot, go->ekin, go->rmsd);
  cago_close(go);
  return 0;
}
