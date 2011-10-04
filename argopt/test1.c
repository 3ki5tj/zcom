#include "argopt.c"

int main(int argc, const char **argv)
{
  argopt_t *ao;
  const char *fn = NULL;
  int help = 0, n = 0, freq = -1, verbose = 0;
  double x = 0;

  ao = argopt_open();
  argopt_regarg(ao, NULL, &fn, "inputfile");
  argopt_regopt(ao, "-h", 0, NULL, &help, "help");
  argopt_regopt(ao, "--verbose", 0, NULL, &verbose, "verbose");
  argopt_regopt(ao, "--freq", 1, "%d", &freq, "frequency of saving files");
  argopt_regopt(ao, "-n", 1, "%d", &n, "integer n");
  argopt_regopt(ao, "-x", 1, "%lf", &x, "double x");
  argopt_parse(ao, argc, argv);
  if (help) argopt_help(ao);
  argopt_close(ao);
  printf("fn %s, n %d, x %g, freq %d\n", fn, n, x, freq);
  return 0;
}
