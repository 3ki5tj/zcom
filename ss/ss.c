#ifndef SS__
#define SS__

#include "ss.h"

#ifndef SSMINSIZ /* to override the block size, define it before inclusion */
#define SSMINSIZ 256 /* at least sizeof(int), but change this value to 1 for debugging */
#endif
#ifndef SSHASHBITS
#define SSHASHBITS 8
#endif
#define SSHASHSIZ   (1 << SSHASHBITS)
#define SSOVERALLOC 1
#define sscalcsize_(n) (((n)/SSMINSIZ + 1) * SSMINSIZ) /* size for n nonblank characters */

struct ssheader {
  size_t size;
  size_t hashval;
  struct ssheader *next;
} ssbase_[SSHASHSIZ] = {{ 0u, 0u, NULL }};

/* we use the string address instead of that of the pointer
 * to struct ssheader to compute the Hash value,
 * because the former is more frequently used in e.g. looking-up
 * */
INLINE size_t sshashval_(const char *p)
{
  size_t val = (size_t) p * 1664525u + 1013904223u;
  return (val >> (sizeof(size_t)*8-SSHASHBITS)) & ((1<<SSHASHBITS)-1);
}

/* return the *previous* header to the one that associates with s
 * first locate the list from the Hash value, then enumerate the linked list.
 * */
INLINE struct ssheader *sslistfind_(const char *s)
{
  struct ssheader *hp, *head;

  if (s == NULL) return NULL;
  head = ssbase_ + sshashval_(s);
  if (head->next == NULL) return NULL; /* uninitialized head node */
  for (hp = head; hp->next != head; hp = hp->next)
    if ((char *)(hp->next + 1) == s) /* match the string address */
      return hp;
  return NULL;
}

/* simply add the entry h at the begining of the list
 * we do not accept a precalculated hash value,
 * since realloc might have changed it
 * */
INLINE struct ssheader *sslistadd_(struct ssheader *h)
{
  struct ssheader *head;

  head = ssbase_ + sshashval_( (char *)(h+1) );
  if (head->next == NULL) /* initialize the base */
    head->next = head;
  h->next = head->next;
  head->next = h;
  return head;
}

/* remove hp->next */
INLINE void sslistremove_(struct ssheader *hp, int f)
{
  struct ssheader *h = hp->next;

  hp->next = h->next;
  if (f) free(h);
}

/* (re)allocate memory for (*php)->next, update list, return the new string
 * n is the number of nonempty characters, obtained e.g. from strlen().
 * create a new header if *php is NULL, in this case, the first character
 * of the string is '\0'
 * */
INLINE char *ssresize_(struct ssheader **php, size_t n, unsigned flags)
{
  struct ssheader *h = NULL, *hn, *hp;
  size_t size;

  if (php == NULL) {
    fprintf(stderr, "ssresize_: php is NULL, n = %u", (unsigned) n);
    return NULL;
  }

  /* we use the following if to assign hp and h, so the order is crucial */
  if ((hp = *php) == NULL || (h = hp->next)->size < n + 1 || !(flags & SSOVERALLOC)) {
    size = sscalcsize_(n);
    if (h == NULL || size != h->size) {
      /* remove h from the list  */
      if (hp != NULL) sslistremove_(hp, 0);
      if ((hn = calloc(sizeof(*h) + size, 1)) == NULL) {
        fprintf(stderr, "ssresize_: no memory for %u\n", (unsigned) size);
        return NULL;
      }
      if (h != NULL) {
        memcpy(hn, h, sizeof(*hn) + (size > h->size ? h->size : size));
        free(h);
      }
      h = hn;

      *php = hp = sslistadd_(h);
      hp->next->size = size;
    }
  }
  return (char *)(hp->next + 1);
}

INLINE void ssmanage_low_(struct ssheader *hp, unsigned opt)
{
  if (opt == SSDELETE)
    sslistremove_(hp, 1);
  else if (opt == SSSHRINK)
    ssresize_(&hp, strlen((char *)(hp->next+1)), 0);
}

/* delete a string, shrink memory, etc ... */
INLINE int ssmanage(char *s, unsigned flags)
{
  struct ssheader *hp, *head;
  unsigned opt = flags & 0xFF;
  size_t i;

  if (flags & SSSINGLE) { /* working on a single string */
    if (s == NULL || (hp = sslistfind_(s)) == NULL) {
      if (s) fprintf(stderr, "ssmanage: unknown address %p (%s)\n",  s, s);
      return -1;
    }
    ssmanage_low_(hp, opt);
  } else {
    for (i = 0; i < SSHASHSIZ; i++)
      for (hp = head = ssbase_+i; hp->next && hp->next != head; hp = hp->next)
        /* we must not operate on h itself, which renders the iterator h invalid */
        ssmanage_low_(hp, opt);
  }
  return 0;
}

/* copy/cat t to *ps
 *
 * If (flags & SSCAT) == 0:
 * copy t to *ps, if ps is not NULL, and return the result
 * if ps or *ps is NULL, we return a string created from t
 *   *ps is set to the same value if ps is not NULL
 * otherwise, we update the record that corresponds to *ps
 *
 * minsize: to request a minimal size for the resulting buffer
 *
 * If flags & SSCAT:
 * append t after *ps. Equivalent to cpy if ps or *ps is NULL.
 * */
INLINE char *sscpycatx(char **ps, const char *t, size_t minsize, unsigned flags)
{
  struct ssheader *hp = NULL;
  size_t size = 0u, sizes = 0u;
  char *s = NULL, *p;

  /* both ps and *ps can be NULL, in which cases we leave hp as NULL */
  if (ps != NULL && (s = *ps) != NULL && (hp = sslistfind_(s)) == NULL) {
    fprintf(stderr, "sscpycatx: unknown address %p (%s)\n", s, s);
    return NULL;
  }
  if (t != NULL)
    while (t[size]) /* compute the length of t */
      size++;
  if (flags & SSCAT) {
    if (s != NULL)  /* s is also NULL, if ps itself is NULL */
      while (s[sizes]) /* compute the length of s */
        sizes++;
    size += sizes;
  }  /* sizes is always 0 in case of copying */
  if (size < minsize)
    size = minsize;
  if ((s = ssresize_(&hp, size, SSOVERALLOC)) == NULL) { /* change size */
    return NULL;
  }
  if (t != NULL)
    for (p = s + sizes; (*p++ = *t++) != '\0'; ) /* copy/cat the string */
      ;
  if (ps != NULL)
    *ps = s;
  return s;
}

/* get a string *ps from file fp
 * *ps can be NULL, in which case memory is allocated
 * *pn is number of characters read (including '\n', but not the terminal null)
 * delim is the '\n' for reading a singe line
 * */
INLINE char *ssfgetx(char **ps, size_t *pn, int delim, FILE *fp)
{
  size_t n, max;
  int c;
  char *s;
  struct ssheader *hp;

  if (ps == NULL || fp == NULL)
    return NULL;
  if ((s = *ps) == NULL) /* allocate an initial buffer if *ps is NULL */
    if ((s = sscpycatx(ps, NULL, 0, 0u)) == NULL)
      return NULL;
  if ((hp = sslistfind_(s)) == NULL) {
    fprintf(stderr, "ssfgetx: unknown address %p (%s)\n", s, s);
    return NULL;
  }
  max = hp->next->size-1;
  for (n = 0; (c = fgetc(fp)) != EOF; ) {
    if (n+1 > max) { /* request space for n+1 nonblank characters */
      if ((*ps = s = ssresize_(&hp, n+1, SSOVERALLOC)) == NULL)
        return NULL;
      max = hp->next->size - 1;
    }
    s[n++] = (char)(unsigned char) c;
    if (c == delim)
      break;
  }
  s[n] = '\0';
  if (pn != NULL)
    *pn = n;
  return (n > 0) ? s : NULL;
}


/* parse `s' into a string array
 * delimiters are removed */
INLINE char **ssparse(char *s, int *pn, const char *delim)
{
  const int capsz = 16;
  int cap, n;
  char **sarr, *p, *q;
  char delim0[8] = "\n\r"; /* default deliminators: new lines */

  if (pn) *pn = 0;
  if (delim == NULL) delim = delim0;

  cap = capsz;
  if ((sarr = calloc(cap, sizeof(sarr[0]))) == NULL) {
    fprintf(stderr, "no memory for sarr\n");
    return NULL;
  }
  for (n = 0, p = s; ; ) { /* n is # of lines */
    for (q = p; *q != '\0'; q++)
      if (strchr(delim, *q))
        break;
    if (q != p) { /* skip an empty line */
      sarr[n++] = p;
      if (n >= cap) { /* expand the array */
        cap += capsz;
        if ((sarr = realloc(sarr, cap * sizeof(sarr[0]))) == NULL) {
          fprintf(stderr, "no memory for sarr, %d\n", cap);
          return NULL;
        }
      }
    }
    if (*q == '\0') break; /* we are done */
    *q = '\0';
    /* search for the next starting point */
    for (p = q + 1; *p && strchr(delim, *p); p++)
      ;
    if (*p == '\0') break;
  }

  if (pn) *pn = n;
  return sarr;
}

/* free the string array, sarr[0] created by ssnew() and sarr created by malloc() */
#define ssarrfree(sarr) { ssdel(sarr[0]); free(sarr); }

#endif

