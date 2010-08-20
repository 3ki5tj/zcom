#ifndef SAFESTRING
#define SAFESTRING 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "ss.h"

#define SSMINBITS  8   /* change this value to 0 during debug */
#define SSMINSIZ   (1u<<SSMINBITS)
#define SSHASHBITS 8
#define SSHASHSIZ  (1<<SSHASHBITS)  

typedef struct tag_ssinfo{
  size_t size;
  size_t hashval;
  struct tag_ssinfo *next;
}ssinfo_t;

static struct tag_ssinfo ssbase_[SSHASHSIZ]={{0u,0u,NULL}};

#define sserror_ printf

/* we accept a string instead of a point to ssinfo_t
 * because the former is more frequently used in e.g. looking-up
 * */
static size_t sshashval_(const char *p)
{
  size_t val=(size_t)p;
  val = val*1664525u+1013904223u;
  return (val >> (sizeof(size_t)*8-SSHASHBITS)) & ((1<<SSHASHBITS)-1);
}

/* calculate memory needed for a string of n nonzero characters */
static size_t sscalcsize_(size_t n, int ah)
{
  n = (n/SSMINSIZ + 1) * SSMINSIZ;
  if (ah) n += sizeof(ssinfo_t);
  return n;
}

/* slow but safe way of verifying the validity of a string
 * by enumerate the whole linked list 
 * return the *previous* item to the requested one */
static ssinfo_t *sslistfind_(char *s)
{
  ssinfo_t *hp;

  if (s == NULL)
    return NULL;
  for (hp = ssbase_ + sshashval_(s); hp->next != ssbase_; hp = hp->next)
    if ((char *)(hp->next+1) == s)
      return hp;
  return NULL;
}

/* simply add h to the begining of the list 
 * we do not accept a precalculated hash value, 
 * since realloc might have changed it
 * */
static ssinfo_t *sslistadd_(ssinfo_t *h)
{
  ssinfo_t *head;

  head = ssbase_ + sshashval_( (char *)(h+1) );
  if (head->next == NULL) /* initialize the base */
    head->next = head;

  h->next = head->next;
  head->next = h;
  return head;
}

/* remove hp->next */
static void sslistremove_(ssinfo_t *hp, int f)
{
  ssinfo_t *h = hp->next;
  hp->next = h->next;
  if (f) free(h);
}

/* (re)allocate memory for (*php)->next, update list, return the new string
 * n is the number of nonempty characters, obtained e.g. from strlen()
 * we create a new header if *php is NULL
 * */
static char *ssresize_(ssinfo_t **php, size_t n, int overalloc)
{
  ssinfo_t *h=NULL, *hp;
  size_t size;

  if (php == NULL)
    sserror_("NULL pointer to resize");
  
  /* we use the following if to assign hp and h, so the order is crucial */
  if ((hp=*php) == NULL || (h = hp->next)->size < (n+1) || !overalloc) {
    size = sscalcsize_(n, 0);
    if (h == NULL || size != h->size) {
      /* since realloc will change the hash value of h
       * we have to remove the old entry first without free() 
       * hp->next will be freed by realloc */
      if (hp != NULL)
        sslistremove_(hp, 0);
      if ((h = realloc(h, sizeof(*h)+size)) == NULL) {
        sserror_("no memory for resizing\n");
        return NULL;
      }
      if (hp == NULL) /* clear the first byte if we start from nothing */
        *(char *)(h+1) = '\0';
      *php = hp = sslistadd_(h);
      hp->next->size = size;
    }
  }
  return (char *)(hp->next+1);
}

static void ssmanage_low_(ssinfo_t *hp, unsigned opt)
{
  if (opt == SSDELETE)
    sslistremove_(hp, 1);
  else if (opt==SSSHRINK)
    ssresize_(&hp, strlen((char *)(hp->next+1)), 0);
  else
    sserror_("unknown manage option");
}

/* delete a string, shrink memory, etc ... */
void ssmanage(char *s, unsigned flags)
{
  ssinfo_t *hp,*head;
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

/* copy t to *ps, if ps is not NULL, and return the result
 * if ps or *ps is NULL, we return a string created from t
 *   *ps is set to the same value if ps is not NULL
 * otherwise, we update the record that corresponds to *ps
 * */
char *sscpyx(char **ps, const char *t)
{
  ssinfo_t *hp=NULL;
  size_t size=0;
  char *s, *p;

  if (ps != NULL && (s=*ps) != NULL && (hp = sslistfind_(s)) == NULL)
    return NULL;
  if (t != NULL)
    while(t[size])
      size++;
  if ((s = ssresize_(&hp, size, 1)) == NULL)
    return NULL;
  if (t == NULL)
    s[0] = '\0';
  else /* copy t to s */
    for (p = s; (*p++ = *t++); )
      ;
  if (ps != NULL)
    *ps = s;
  return s;
}

char *sscatx(char **ps, const char *t)
{
  ssinfo_t *hp;
  size_t size=0;
  char *s, *p;

  if (ps == NULL || (s=*ps) == NULL || (hp = sslistfind_(s)) == NULL || t == NULL)
    return NULL;
  while (t[size])
    size++;
  for (p = s; *p; p++) ; /* move p to the end of the string */
  size += p-s;
  if ((s = ssresize_(&hp, size, 1)) == NULL)
    return NULL;
  if (s != *ps) /* move to the end of s */
    for (p = *ps = s; *p; p++)
      ;
  for (; (*p++ = *t++); )
    ;
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
  ssinfo_t *hp;

  if (ps == NULL || fp == NULL)
    return NULL;
  if ((s=*ps) == NULL) /* allocate an initial buffer if *ps is NULL */
    if ((s = sscpyx(ps,NULL)) == NULL)
      return NULL;
  if ((hp = sslistfind_(s)) == NULL)
    return NULL;
  max = hp->next->size-1;
  for (n = 0; (c = fgetc(fp)) != EOF; ) {
    if (n+1 > max) { /* request space for n+1 nonblank characters */
      if ((*ps = s = ssresize_(&hp, n+1, 1)) == NULL)
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

