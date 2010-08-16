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
#ifdef SSDBG_
static void ssdump_(ssinfo_t *hp)
{
  fflush(stdout);
  printf("header: %p, prev: %p, next: %p, base: %p, size: %8u, string: %p\n",
      hp->next, hp, hp->next->next, ssbase_, hp->next->size, (char *)(hp->next+1));
  fflush(stdout);
}
static void ssdumpall_(const char *fmt, ...)
{
  ssinfo_t *hp;
  int i;
  va_list args;
  
  if(fmt){
    fflush(stdout);
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
  fflush(stdout);
  for(hp=ssbase_, i=0; hp->next != ssbase_; hp=hp->next, i++){
    printf("i=%4d:", i);
    ssdump_(hp);
  }
}
#endif

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
#ifdef SSDBG_
    if(hp != NULL){
        printf("this is prechecking, hp=%p size=%u\n", hp, hp->size); fflush(stdout);
        printf("hp->next=%p, size=%u\n",
            hp->next, hp->next->size); fflush(stdout);
        printf("hp->next->next=%p, size=%u\n",
            hp->next->next, hp->next->next->size); fflush(stdout);
    }
#endif
    
#ifdef SSDBG_
    if(hp != NULL){
      printf("h=%p ", h); fflush(stdout);
      if(h==NULL){
        ssdumpall_("hp=%p\n", hp);
        exit(1);
      }
      printf("current size=%u ", h->size); fflush(stdout);
    }
    printf("request %u after %u, overalloc=%d\n", size, sscalcsize_(size,0), overalloc);
#endif
    size = sscalcsize_(size, 0);
    if(h == NULL || size != h->size){
      if((h=realloc(h, sizeof(*hp)+size)) == NULL){
        sserror_("no memory for resizing\n");
        return NULL;
      }
      if(hp != NULL){ /* no need to change h's position in the list */
#ifdef SSDBG_
        printf("this is re-sizing, hp=%p size=%u\n", hp, hp->size); fflush(stdout);
        printf("h=%p, size=%u\n", h, h->size); fflush(stdout);
#endif
        hp->next = h;
      }else{ /* absorb the new header */
#ifdef SSDBG_
        printf("this is a new object\n");
#endif
        *php=hp=sslistadd_(h);
      }
      hp->next->size = size;
#ifdef SSDBG_
      ssdumpall_(NULL); 
#endif
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
#ifdef SSDBG_    
    ssdumpall_("after ssmanage single mode\n");
#endif
  }else{
    for(hp=ssbase_; hp->next != ssbase_; hp=hp->next){
      /* we must not operate on h itself, which renders the iterator h invalid */
#ifdef SSDBG__
      printf("hahaha hp=%p, hp->next=%p, base=%p ", hp, hp->next, ssbase_); fflush(stdout);
      printf("hp->next->next=%p, hp->next:s=%u\n", 
          hp->next->next, strlen((char *)(hp->next+1)) ); getchar();
#endif
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
#ifdef SSDBG_
  if(t){
    printf("sscpyx: size of t=%u, strlen=%u\n", size, strlen(t));
    if(strlen(t) != size){
      fprintf(stderr, "Fatal error\n");
      exit(1);
    }
  }
  if(ps && *ps){
    printf("sscpyx: cap=%u, size=%u\n", hp->next->size, strlen(*ps));
  }
#endif
  if((s=ssresize_(&hp, size, 1)) == NULL)
    return NULL;
#ifdef SSDBG_
  printf("sscpyx: after resizing, cap=%u\n", hp->next->size);
#endif
  if(t == NULL)
    s[0] = '\0';
  else /* copy t to s */
    for(p=s; (*p++ = *t++); )
      ;
#ifdef SSDBG_
  if(t){
    printf("sscpyx: p-s=%u\ns=%s\n", p-s, s); getchar();
  }
#endif
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
#ifdef SSDBG_
  printf("sscatx: size of t=%u\n", size);
#endif
  for(p=s; *p; p++) ; /* move p to the end of the string */
  size += p-s;
#ifdef SSDBG_
  printf("sscatx: size of s=%u, p-s=%u, cap=%u\ns=%s\nt=%s\n",
      size, p-s, hp->next->size, s, t);
#endif
  if((s=ssresize_(&hp, size, 1)) == NULL)
    return NULL;
#ifdef SSDBG_
  printf("sscatx: after resizing hp=%p ", hp); fflush(stdout); 
  printf("h=%p ", hp->next); fflush(stdout);
  printf("s=%p, %p ", hp->next+1, s); fflush(stdout);
  printf("*s=%c\n", *s); fflush(stdout);
  printf("sscatx: after resizing s=%s\naddr: %p\nprev: %p\ncap=%u\n",
      s, s, *ps, hp->next->size); getchar();
#endif
  if(s != *ps) /* move to the end of s */
    for(p=*ps=s; *p; p++)
      ;
#ifdef SSDBG_
  printf("sscatx: p-s=%u\n", p-s);
#endif
  for(; (*p++ = *t++); )
    ;
#ifdef SSDBG_
  printf("sscatx: p-s=%u\ns=%s\n", p-s, s); getchar();
#endif
  return s;
}

#endif

