#include "argopt.h"

int main(int argc, char **argv)
{
  argopt_t *ao;
  const char *fn = NULL;
  char *str = NULL;
  int n = 0, freq = -1, verbose = 0;
  double x = 1e-31;
  real y = 0;
  int nset, xset, yset, fset;

  ao = argopt_open(0);

  argopt_add(ao, NULL, "!", &fn, "input file");
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
  argopt_add(ao, "-s", NULL, &str, "a string (allowing , as the separator)");
  argopt_addversion(ao, "--version");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  nset = argopt_isset(ao, n);
  xset = argopt_isset(ao, x);
  yset = argopt_isset(ao, y);
  fset = argopt_isset(ao, freq);
  printf("fn %s, verbose %d, n[%d] %d, x[%d] %g, y[%d] %g, freq[%d] %d, str ``%s''\n",
      fn, verbose, nset, n, xset, x, yset, y, fset, freq, str);
  argopt_close(ao);
  return 0;
}
