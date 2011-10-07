#include "lj.c"
#include "include/av.h"

int N = 256;
real rho = 0.8f;
real rc = 2.5f;
real tp = 1.0f;
real mddt = 0.002f;
real thermdt = 0.01f;
av_t avU, avK, avp;

int main(void)
{
  lj_t *lj;
  int t, nsteps = 20000;
  
  lj = lj_open(N, 3, rho, rc, tp);
  for (t = 0; t < nsteps; t++) {
    lj_vv(lj, mddt);
    lj_vrescale(lj, thermdt);
    if (t > nsteps/2) {
      av_add(&avU, lj->epot);
      av_add(&avK, lj->ekin);
      av_add(&avp, lj->p);
    }
  }
  printf("U/N = %6.3f, K/N = %6.3f, p = %6.3f\n", av_getave(&avU)/N, av_getave(&avK)/N, av_getave(&avp));  
  lj_close(lj);
  return 0;
}
