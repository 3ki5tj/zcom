#ifndef SAFESTRING
#define SAFESTRING 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "ss.h"

#define SSMINSIZ_BITS 8
#define SSMINSIZ      (1u<<SSMINSIZ_BITS)

typedef struct tag_ssinfo{
  size_t size;
  struct tag_ssinfo *next;
}ssinfo_t;

static struct tag_ssinfo ssbase_[1]={{0u,NULL}};

#define sserror_ printf

/* calculate memory required to accommendate a string of n nonzero characters  */
static size_t sscalcsize_(size_t n, int ah)
{
  n = (n/SSMINSIZ + 1) * SSMINSIZ;
  if(ah) n+=sizeof(ssinfo_t);
  return n;
}

/* slow but safe way of verifying the validity of a string
 * by enumerate the whole linked list 
 * return the *previous* item to the requested one */
static ssinfo_t *sslistfind_(char *s)
{
  ssinfo_t *hp;

  for(hp=ssbase_; hp->next != ssbase_; hp=hp->next)
    if((char *)(hp->next+1) == s)
      return hp;
  return NULL;
}

/* simply add h to the begining of the list */
static ssinfo_t *sslistadd_(ssinfo_t *h)
{
  if(ssbase_->next == NULL) /* initialize the base */
    ssbase_->next = ssbase_;

  h->next = ssbase_->next;
  ssbase_->next = h;
  return ssbase_;
}

/* remove hp->next */
static void sslistremove_(ssinfo_t *hp)
{
  ssinfo_t *h=hp->next;
  hp->next=h->next;
  free(h);
}

/* (re)allocate memory for (*php)->next, update list, return the new string
 * we create a new header if *php is NULL
 * the allocate string will be no less than the requested size
 * */
static char *ssresize_(ssinfo_t **php, size_t size, int overalloc)
{
  ssinfo_t *h=NULL, *hp;

  if(php == NULL){
    sserror_("NULL pointer to resize");
    exit(1);
  }
  /* we use the following if to assign hp and h, so the order is crucial */
  if((hp=*php) == NULL || (h=hp->next)->size < size || !overalloc){
    
    size = sscalcsize_(size, 0);
    if(h == NULL || size != h->size){
      if((h=realloc(h, sizeof(*hp)+size)) == NULL){
        sserror_("no memory for resizing\n");
        return NULL;
      }
      if(hp != NULL){ /* no need to change h's position in the list */
        hp->next = h;
      }else{ /* absorb the new header */
        *php=hp=sslistadd_(h);
      }
      hp->next->size = size;
    }
  }
  return (char *)(hp->next+1);
}

static void ssmanage_low_(ssinfo_t *hp, unsigned opt)
{
  if(opt == SSDELETE)
    sslistremove_(hp);
  else if(opt==SSSHRINK)
    ssresize_(&hp, strlen((char *)(hp->next+1)), 0);
  else
    sserror_("unknown manage option");
}

/* delete a string, shrink memory, etc ... */
void ssmanage(char *s, unsigned flags)
{
  ssinfo_t *hp;
  unsigned opt = flags & 0xFF;

  if(flags & SSSINGLE){
    if(s == NULL || (hp=sslistfind_(s)) == NULL)
      return;
    ssmanage_low_(hp, opt);
  }else{
    for(hp=ssbase_; hp->next != ssbase_; hp=hp->next){
      /* we must not operate on h itself, which renders the iterator h invalid */
      ssmanage_low_(hp, opt);
    }
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

  if(ps != NULL && (s=*ps) != NULL && (hp=sslistfind_(s)) == NULL)
    return NULL;
  if(t != NULL)
    while(t[size])
      size++;
  if((s=ssresize_(&hp, size, 1)) == NULL)
    return NULL;
  if(t == NULL)
    s[0] = '\0';
  else /* copy t to s */
    for(p=s; (*p++ = *t++); )
      ;
  if(ps != NULL)
    *ps = s;
  return s;
}

char *sscatx(char **ps, const char *t)
{
  ssinfo_t *hp;
  size_t size=0;
  char *s, *p;

  if(ps == NULL || (s=*ps) == NULL || (hp=sslistfind_(s)) == NULL || t == NULL)
    return NULL;
  while(t[size])
    size++;
  for(p=s; *p; p++) ; /* move p to the end of the string */
  size += p-s;
  if((s=ssresize_(&hp, size, 1)) == NULL)
    return NULL;
  if(s != *ps) /* move to the end of s */
    for(p=*ps=s; *p; p++)
      ;
  for(; (*p++ = *t++); )
    ;
  return s;
}

#endif

