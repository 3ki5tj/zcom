#include "argopt.c"

int main(int argc, char **argv)
{
  argopt_t *ao;
  const char *fn = NULL;
  int n = 0, freq = -1, verbose = 0;
  double x = 1e-31;
  real y = 0;

  ao = argopt_open(0);


  argopt_add(ao, NULL, "!", &fn, "inputfile");
/*
  the ``NULL'' means argument, equivalent to

    argopt_regarg(ao, "!", &fn, "inputfile");

  "!" means a must
*/

  /* argopt_add() == argopt_regopt() */
  argopt_add(ao, "--verbose", "%b", &verbose, "verbose");
  argopt_add(ao, "--freq", "%d", &freq, "frequency of saving files");
  argopt_add(ao, "-n", "!%d", &n, "an integer n");
  argopt_add(ao, "-x", "%lf", &x, "a double x");
  argopt_add(ao, "-y", "%r", &y, "a real y");
  argopt_addversion(ao, "--version");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  printf("fn %s, verbose %d, n[%d] %d, x[%d] %g, y[%d] %g, freq[%d] %d\n",
      fn, verbose, argopt_set(ao, n),  n, argopt_set(ao, x), x, argopt_set(ao, y), y, 
      argopt_set(ao, freq), freq);
  argopt_close(ao);
  return 0;
}
