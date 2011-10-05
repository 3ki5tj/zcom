#include "argopt.c"

int main(int argc, const char **argv)
{
  argopt_t *ao;
  const char *fn = NULL;
  int n = 0, freq = -1, verbose = 0;
  double x = 0;
  real y = 0;
  opt_t *on, *ox, *oy, *ofreq;

  ao = argopt_open(0, "Test program", "James Bond");
  argopt_regarg(ao, 1, NULL, &fn, "inputfile");
  argopt_regopt_help(ao, "-h");
  argopt_regopt_version(ao, "--version");
  argopt_regopt(ao, "--verbose", 0, NULL, &verbose, "verbose");
  ofreq = argopt_regopt(ao, "--freq", 1, "%d", &freq, "frequency of saving files");
  on = argopt_regopt(ao, "-n", 1, "%d", &n, "integer n");
  ox =  argopt_regopt(ao, "-x", 1, "%lf", &x, "double x");
  oy = argopt_regopt(ao, "-y", 1, "%r", &y, "\'real\' y");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
  printf("fn %s, n[%d] %d, x[%d] %g, y[%d] %g, freq[%d] %d\n",
      fn, on->set, n, ox->set, x, oy->set, y, ofreq->set, freq);
  return 0;
}
