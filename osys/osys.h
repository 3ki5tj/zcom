#include "util.c"
#include "ss.c"

#ifndef OSYS_H__
#define OSYS_H__

#define OSYS_VERBOSE 0x1

/* run command and capture the output
 * `cmd' is the command to run
 * `*nc' returns the number of characters read */
INLINE char *sysrun(const char *cmd, size_t *nc, unsigned flags)
{
  FILE *fp;
  char fntmp[] = "h1Qi7G0c.TmP"; /* output file */
  char *ncmd, *output = NULL;
  int i;

  if (nc) *nc = 0;

  /* construct the command */
  ncmd = ssdup(cmd);
  sscat(ncmd, " > ");
  sscat(ncmd, fntmp);

  /* run it */
  if ((i = system(ncmd)) != 0)
    if (flags & OSYS_VERBOSE)
      fprintf(stderr, "command \"%s\" failed\n", ncmd);
  ssdel(ncmd);
  if (i != 0) goto EXIT; /* output should be NULL in this case */

  xfopen(fp, fntmp, "r", goto EXIT);
  ssfgetall(output, nc, fp);
  fclose(fp);
EXIT:
  remove(fntmp);
  return output;
}


/* Get a list of file names satisfying a pattern
 * `pat' is the file name pattern, such as '*.c'
 * `*n' returns the number of matches
 * `lscmd' is the list command
 * the output is a file list (string array), to free it:
 *    ssdel(fnls[0]); free(fnls); */
INLINE char **fnglob(const char *pat, int *pn, const char *lscmd, unsigned flags)
{
  size_t nc = 0;
  char *cmd, *output;

  if (pn) *pn = 0;

  /* construct a command first */
  cmd = ssnew(256);
  if (lscmd == NULL) /* assuming Linux */
    lscmd = "ls --color=never ";
  sscpy(cmd, lscmd);
  sscat(cmd, pat);

  /* run the command and get the output */
  output = sysrun(cmd, &nc, flags);
  ssdel(cmd);
  if (output == NULL || nc == 0) return NULL;

  /* parse the output into file names
   * don't ssdel(output), for it's used in the output */
  return ssparse(output, pn, NULL);
}


#endif

