#include "argopt.c"

int main(int argc, char **argv)
{
  argopt_t *ao;
  const char *fn = NULL;
  int n = 0, freq = -1, verbose;
  double x = 0;
  real y = 0;

  ao = argopt_open(0);
  argopt_regarg(ao, "!", &fn, "inputfile");
  argopt_regopt_help(ao, "-h");
  argopt_regopt_version(ao, "--version");
  argopt_regopt(ao, "--verbose", "%b", &verbose, "verbose");
  argopt_regopt(ao, "--freq", "%d", &freq, "frequency of saving files");
  argopt_regopt(ao, "-n", "!%d", &n, "an integer n");
  argopt_regopt(ao, "-x", "%lf", &x, "a double x");
  argopt_regopt(ao, "-y", "%r", &y, "a real y");
  argopt_parse(ao, argc, argv);
  printf("fn %s, verbose %d, n[%d] %d, x[%d] %g, y[%d] %g, freq[%d] %d\n",
      fn, verbose, argopt_set(ao, n),  n, argopt_set(ao, x), x, argopt_set(ao, y), y, 
      argopt_set(ao, freq), freq);
  argopt_close(ao);
  return 0;
}