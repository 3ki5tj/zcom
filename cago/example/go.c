#define ZCOM_PICK
#define ZCOM_CAGO
#define ZCOM_ARGOPT
#define ZCOM_HIST
#define ZCOM_AV
#include "zcom.h"

/* Ca-Go model parameters */
real kb = 200.f;
real ka = 40.f;
real kd1 = 1.f;
real kd3 = 0.5f;
real nbe = 1.f;
real nbc = 4.f;
real rcc = 4.5f;

int nsteps = 1000000000;
real mddt = 0.002f;
real tp = 1.12f;

real thermdt = 0.1f; /* thermostat dt */
int method = 0;


/* for SH3 domain (1KIK), if rcc = 4.5 --> then Tc = 1.12 
 * the double-peak structure can be seen from the potential-energy
 * and rmsd distributions */
char *fnpdb = "pdb/1KIK.pdb";
real umin = 0.f;
real umax = 140.f;
real udel = 2.0f;
/* parameters for the potential energy histogram */
real uhmin = -300.f;
real uhmax = 300.f;
real uhdel = 0.1f;
int nevery = 10000;
int nreport = 500000;
int nequil = 10000;
const char *fnephis = "epot.his";
const char *fnrmsdhis = "rmsd.his";
const char *fnconthis = "cont.his";
const char *fnpos = "go.pos";

/* handle input arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  argopt_add(ao, "-f", NULL, &fnpdb,    "PDB file");
  argopt_add(ao, "-T", "%r", &tp,       "default temperature");
  argopt_add(ao, "-1", "%d", &nsteps,   "number of simulation steps");
  argopt_add(ao, "-d", "%r", &mddt,     "time step for molecular dynamics");
  argopt_add(ao, "-q", "%r", &thermdt,  "time step for mc-vrescaling thermostat");
  argopt_add(ao, "-c", "%r", &rcc,      "cutoff distance for defining contacts");
  argopt_add(ao, "-m", "%d", &method,   "0: mc samp; 1: vrescale");  
  argopt_add(ao, "--umin", "%r", &umin,     "minimal potential energy");
  argopt_add(ao, "--umax", "%r", &umax,     "maximal potential energy");
  argopt_add(ao, "--udel", "%r", &udel,     "potential energy interval");
  argopt_add(ao, "--uhmin", "%r", &uhmin,   "minimal potential energy for histogram");
  argopt_add(ao, "--uhmax", "%r", &uhmax,   "maximal potential energy for histogram");
  argopt_add(ao, "--uhdel", "%r", &uhdel,   "potential energy interval for histogram");
  argopt_add(ao, "--every", "%d", &nevery,   "print message every this # of steps");
  argopt_add(ao, "--report", "%d", &nreport, "save data every this # of steps");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);

  argopt_dump(ao);
  argopt_close(ao);
}

static void domd(void)
{
  cago_t *w;
  int it, ncont;
  real rmsd;
  av_t avep[1] = {{0}}, avrmsd[1] = {{0}}, avcont[1] = {{0}};
  hist_t *hsep, *hsrmsd, *hscont;

  w = cago_open(fnpdb, kb, ka, kd1, kd3, nbe, nbc, rcc);
  cago_initmd(w, 0.1, 0.0);
  printf("epot %g, epotref %g\n", w->epot, w->epotref);

  hsep = hs_open(1, uhmin, uhmax, uhdel);
  hsrmsd = hs_open(1, 0, 100.0, 0.1);
  hscont = hs_open(1, 0, w->ncont + 1, 1);

  for (it = 1; it <= nsteps; it++) {
    cago_vv(w, 1.0f, mddt);
    cago_rmcom(w, w->x, w->v);
    cago_vrescale(w, tp, thermdt);
    rmsd = cago_rmsd(w, w->x, NULL); /* compute rmsd */
    ncont = cago_countcontact(w, w->x, 1.2, NULL, NULL);
    
    if (it >= nequil) {
      av_add(avep, w->epot);
      av_add(avrmsd, rmsd);
      av_add(avcont, ncont); 
      hs_add1ez(hsep, w->epot, HIST_VERBOSE); /* add to histogram */
      hs_add1ez(hsrmsd, rmsd, HIST_VERBOSE);
      hs_add1ez(hscont, ncont, HIST_VERBOSE);

      if (it % nevery == 0) {
        printf("t %d, T %g, ep %g/%g, rmsd %g/%g Q %d/%d=%g(%g)\n", it, w->tkin,
          w->epot, av_getave(avep), rmsd, av_getave(avrmsd),
          ncont, w->ncont, 1.0*ncont/w->ncont, av_getave(avcont)/w->ncont);
      }
      if (it % nreport == 0 || it == nsteps) {
        printf("saving %s, %s, %s\n", fnephis, fnrmsdhis, fnpos);
        cago_writepos(w, w->x, w->v, fnpos);
        hs_save(hsep, fnephis, HIST_ADDAHALF);
        hs_save(hsrmsd, fnrmsdhis, HIST_ADDAHALF);
        hs_save(hscont, fnconthis, 0);
      }
    }
  }
  cago_close(w);
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

