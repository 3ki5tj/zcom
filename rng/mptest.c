#include "rng.h"
#include <time.h>
#include <omp.h>

int main(void)
{
#pragma omp parallel
  {
    int i = omp_get_thread_num();

    mtscramble(i + 1 + time(NULL));
    printf("thread %d, mr_ %p, %u mrstock_ %p, %u\n",
        i, (void *) mr_, mr_->arr[0], (void *) &mrstock_, mrstock_.arr[0]);
  }
  return 0;
}
