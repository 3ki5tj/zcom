#define ZCOM_PICK
#define ZCOM_CAGO
#define ZCOM_ARGOPT
#define ZCOM_HIST
#include "zcom.h"

/* Ca-Go model parameters */
real kb  = (real) 200.0;
real ka  = (real) 40.0;
real kd1 = (real) 1.0;
real kd3 = (real) 0.5;
real nbe = (real) 1.0;
real nbc = (real) 4.0;
real rcc = (real) 4.5;
int fromrand = 0;
/* in counting formed contacts in a configuration
 * we consider a contact is formed if the pair distance is shorter 
 * than the `ncgam' * the distance in the reference structure */
double ncgam = 1.2;

int nsteps = 1000000000;
real mddt = 0.002f;
real tp = 1.0f;
real thermdt = 0.1f; /* thermostat dt */

/* for SH3 domain (1KIK), if rcc = 4.5 --> then Tc = 1.12 
 * the double-peak structure can be seen from the potential-energy
 * and rmsd distributions */
char *fnpdb = "pdb/1KIK.pdb";
/* parameters for the potential energy histogram */
int nevery = 10000;
int nreport = 500000;
int nequil = 10000;

const char *fnephis = "ep.his";
const char *fnrmsdhis = "rmsd.his";
const char *fnconthis = "nc.his";
const char *fnpos = "go.pos";

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  argopt_add(ao, "-f", NULL,  &fnpdb,     "PDB file");
  argopt_add(ao, "-Z", "%b",  &fromrand,  "start from random configuration");
  argopt_add(ao, "-G", "%lf", &ncgam,     "distance-scaling factor for counting formed contacts");
  argopt_add(ao, "-T", "%r",  &tp,        "default temperature");
  argopt_add(ao, "-1", "%d",  &nsteps,    "number of simulation steps");
  argopt_add(ao, "-d", "%r",  &mddt,      "time step for molecular dynamics");
  argopt_add(ao, "-q", "%r",  &thermdt,   "time step for mc-vrescaling thermostat");
  argopt_add(ao, "-c", "%r",  &rcc,       "cutoff distance for defining contacts");
  argopt_add(ao, "--every",  "%d", &nevery,  "print messages every this number of steps");
  argopt_add(ao, "--report", "%d", &nreport, "save data every this number of steps");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);

  argopt_dump(ao);
  argopt_close(ao);
}

static void domd(void)
{
  cago_t *go;
  int it, nc;
  real rmsd;
  hist_t *hsep, *hsrmsd, *hscont;

  go = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc);
  cago_initmd(go, fromrand ? -0.1 : 0.1, 0.0);
  printf("epot %.2f (ref %.2f)\n", go->epot, go->epotref);

  hsep = hs_open1(go->epotref*2, fabs(go->epotref)*5, 0.1);
  hsrmsd = hs_open1(0, 100.0, 0.1);
  hscont = hs_open1(0, go->ncont + 1, 1);

  for (it = 1; it <= nsteps; it++) {
    cago_vv(go, 1.0f, mddt);
    cago_rmcom(go, go->x, go->v);
    cago_vrescale(go, tp, thermdt);
    rmsd = cago_rmsd(go, go->x, NULL); /* compute rmsd */
    nc = cago_ncontacts(go, go->x, ncgam, NULL, NULL);
    
    if (it >= nequil) {
      hs_add1ez(hsep, go->epot, 0); /* add to histogram */
      hs_add1ez(hsrmsd, rmsd, HIST_VERBOSE);
      hs_add1ez(hscont, nc, HIST_VERBOSE);

      if (it % nevery == 0) {
        printf("t %d, T %.2f, ep %.2f/%.2f, rmsd %.2f/%.2f Q %d/%d=%.2f(%.2f)\n", it, go->tkin,
          go->epot, hs_getave(hsep, 0, NULL, NULL),
          rmsd, hs_getave(hsrmsd, 0, NULL, NULL),
          nc, go->ncont, 1.0*nc/go->ncont,
          hs_getave(hscont, 0, NULL, NULL)/go->ncont);
      }
      if (it % nreport == 0 || it == nsteps) {
        printf("saving %s, %s, %s, %s\n", fnephis, fnrmsdhis, fnconthis, fnpos);
        cago_writepos(go, go->x, go->v, fnpos);
        hs_save(hsep, fnephis, HIST_ADDAHALF);
        hs_save(hsrmsd, fnrmsdhis, HIST_ADDAHALF);
        hs_save(hscont, fnconthis, 0);
      }
    }
  }
  cago_close(go);
  hs_close(hsep);
  hs_close(hsrmsd);
  hs_close(hscont);
}

int main(int argc, char **argv)
{
  doargs(argc, argv);
  domd();
  return 0;
}

