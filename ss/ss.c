#ifndef SS__
#define SS__

#include "ss.h"

#ifndef SSMINSIZ /* to override the block size, define it before inclusion */
#define SSMINSIZ 256 /* change this value to 1 for debugging */
#endif
#ifndef SSHASHBITS
#define SSHASHBITS 8
#endif
#define SSHASHSIZ  (1<<SSHASHBITS)  
#define SSOVERALLOC 1
#define sscalcsize_(n) (((n)/SSMINSIZ + 1) * SSMINSIZ) /* size for n nonblank characters */

#ifdef ZCOM_ERROR
#define sserror_ fatal
#else
/* print an error message and quit */
static void sserror_(char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  fprintf(stderr, "fatal error: ");
  fprintf(stderr, fmt, args);
  va_end(args);
  exit(1);
}
#endif

struct ssheader{
  size_t size;
  size_t hashval;
  struct ssheader *next;
} ssbase_[SSHASHSIZ] = {{0u, 0u, NULL}};

/* we use the string address instead of that of the pointer 
 * to struct ssheader to compute the Hash value,
 * because the former is more frequently used in e.g. looking-up
 * */
static size_t sshashval_(const char *p)
{
  size_t val = (size_t)p * 1664525u + 1013904223u;
  return (val >> (sizeof(size_t)*8-SSHASHBITS)) & ((1<<SSHASHBITS)-1);
}

/* 
 * return the *previous* header to the one that associates with s
 * first locate the list from the Hash value, then enumerate the linked list.
 * */
static struct ssheader *sslistfind_(char *s)
{
  struct ssheader *hp;

  if (s == NULL)
    return NULL;
  for (hp = ssbase_ + sshashval_(s); hp->next != ssbase_; hp = hp->next)
    if ((char *)(hp->next + 1) == s)
      return hp;
  return NULL;
}

/* 
 * simply add the entry h at the begining of the list 
 * we do not accept a precalculated hash value, 
 * since realloc might have changed it
 * */
static struct ssheader *sslistadd_(struct ssheader *h)
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
static void sslistremove_(struct ssheader *hp, int f)
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
static char *ssresize_(struct ssheader **php, size_t n, unsigned flags)
{
  struct ssheader *h=NULL, *hp;
  size_t size;

  if (php == NULL)
    sserror_("ssresize_: NULL pointer to resize");
  
  /* we use the following if to assign hp and h, so the order is crucial */
  if ((hp=*php) == NULL || (h = hp->next)->size < n + 1 || !(flags & SSOVERALLOC)) {
    size = sscalcsize_(n);
    if (h == NULL || size != h->size) {
      /* since realloc will change the hash value of h
       * we have to remove the old entry first without free() 
       * hp->next will be freed by realloc */
      if (hp != NULL)
        sslistremove_(hp, 0);
      if ((h = realloc(h, sizeof(*h)+size)) == NULL) {
        sserror_("ssresize_: no memory for %u\n", size);
        return NULL;
      }
      if (hp == NULL) /* clear the first byte if we start from nothing */
        *(char *)(h + 1) = '\0';  /* h + 1 is the beginning of the string */
      *php = hp = sslistadd_(h);
      hp->next->size = size;
    }
  }
  return (char *)(hp->next + 1);
}

static void ssmanage_low_(struct ssheader *hp, unsigned opt)
{
  if (opt == SSDELETE)
    sslistremove_(hp, 1);
  else if (opt == SSSHRINK)
    ssresize_(&hp, strlen((char *)(hp->next+1)), 0);
  else
    sserror_("unknown manage option");
}

/* delete a string, shrink memory, etc ... */
void ssmanage(char *s, unsigned flags)
{
  struct ssheader *hp,*head;
  unsigned opt = flags & 0xFF;
  size_t i;

  if (flags & SSSINGLE) {
    if (s == NULL || (hp = sslistfind_(s)) == NULL)
      return;
    ssmanage_low_(hp, opt);
  } else {
    for (i=0; i<SSHASHSIZ; i++)
      for (hp = head = ssbase_+i; hp->next && hp->next != head; hp = hp->next)
        /* we must not operate on h itself, which renders the iterator h invalid */
        ssmanage_low_(hp, opt);
  }
}

/* 
 * copy/cat t to *ps
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
char *sscpycatx(char **ps, const char *t, size_t minsize, unsigned flags)
{
  struct ssheader *hp=NULL;
  size_t size=0u, sizes=0u;
  char *s=NULL, *p;

  /* both ps and *ps can be NULL, in which cases we leave hp as NULL */
  if (ps != NULL && (s=*ps) != NULL && (hp = sslistfind_(s)) == NULL) {
    fprintf(stderr, "sscpycatx: string is not previously registered!\n");
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
    for (p = s + sizes; (*p++ = *t++); ) /* copy/cat the string */
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
char *ssfgetx(char **ps, size_t *pn, int delim, FILE *fp)
{
  size_t n, max;
  int c;
  char *s;
  struct ssheader *hp;

  if (ps == NULL || fp == NULL)
    return NULL;
  if ((s=*ps) == NULL) /* allocate an initial buffer if *ps is NULL */
    if ((s = sscpycatx(ps, NULL, 0, 0u)) == NULL)
      return NULL;
  if ((hp = sslistfind_(s)) == NULL)
    return NULL;
  max = hp->next->size-1;
  for (n = 0; (c = fgetc(fp)) != EOF; ) {
    if (n+1 > max) { /* request space for n+1 nonblank characters */
      if ((*ps = s = ssresize_(&hp, n+1, SSOVERALLOC)) == NULL)
        return NULL;
      max = hp->next->size - 1;
    }
    s[n++] = (char)(unsigned char)c;
    if (c == delim) 
      break;
  }
  s[n] = '\0';
  if (pn != NULL)
    *pn = n;
  return (n > 0) ? s : NULL;
}
#endif
