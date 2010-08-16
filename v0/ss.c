#ifndef SAFESTRING
#define SAFESTRING 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ss.h"

#define SSMINSIZ_BITS 8
#define SSMINSIZ      (1u<<SSMINSIZ_BITS)

typedef struct tag_ssinfo{
  size_t size;
  struct tag_ssinfo *prev, *next;
  size_t mark;
}ssinfo_t;

#define SSMARK (size_t)0x1f3d5b79 /* signature mark */

static struct tag_ssinfo ssbase_[1]={{0u,NULL,NULL,SSMARK}};

#define sserror_ printf

/* calculate memory required to accommendate a string of n nonzero characters  */
static size_t sscalcsize_(size_t n, int hdr)
{
  n = (n/SSMINSIZ + 1) * SSMINSIZ;
  if(hdr) n+=sizeof(ssinfo_t);
  return n;
}

static ssinfo_t *sslistfind_(char *s)
{
  ssinfo_t *hdr;

  if(s == NULL) return NULL;
  hdr = (ssinfo_t *)s - 1;

  /* trying to access hdr->mark may lead to an error,
   * but this is a program error anyway */
  if(hdr->mark != SSMARK){
    sserror_("trying to use an alien string\n");
    return NULL;
  }
  return hdr;
}

static void sslistadd_(ssinfo_t *hdr)
{
  if(ssbase_->prev == NULL) /* initialize the base */
    ssbase_->next = ssbase_->prev = ssbase_;
  
  hdr->next = ssbase_->next;
  ssbase_->next->prev = hdr;

  hdr->prev = ssbase_;
  ssbase_->next = hdr;
}

static void sslistremove_(ssinfo_t *hdr)
{
  hdr->prev->next = hdr->next;
  hdr->next->prev = hdr->prev;
  free(hdr);
}

/* (re)allocate memory for *phdr, update list, return the new string
 * we create a new header if *phdr is NULL
 * the allocate string will be no less than the requested size
 * */
static char *ssresize_(ssinfo_t **phdr, size_t size, int overalloc)
{
  ssinfo_t *hdr;

  if(phdr == NULL){
    sserror_("pass NULL pointer\n");
    return NULL;
  }
  if((hdr=*phdr) == NULL || !overalloc || size > hdr->size){
    size = sscalcsize_(size, 0);
    if(hdr == NULL || size != hdr->size){
      if((hdr=realloc(hdr, sizeof(ssinfo_t)+size)) == NULL){
        sserror_("no memory for resizing\n");
        return NULL;
      }
      if(*phdr != NULL)
        hdr->next->prev = hdr->prev->next = hdr;
      else{ /* absorb the new string */
        sslistadd_(hdr);
      }
      *phdr = hdr;
      hdr->size = size;
      hdr->mark = SSMARK;
    }
  }
  return (char *)(hdr+1);
}

#ifdef SSDBG_
static void ssdump_(ssinfo_t *hdr)
{
  printf("header: %p, prev: %p, next: %p, base: %p, size: %8u, string: %p\n", 
      hdr, hdr->prev, hdr->next, ssbase_, hdr->size, (char *)hdr+sizeof(*hdr));
}
#endif

/* manually delete a string */
void ssdelete(char *s)
{
  ssinfo_t *hdr;

  if(s == NULL || (hdr=sslistfind_(s)) == NULL)
    return;
  sslistremove_(hdr);
} 

/* 0: free everything; 1: shrink */
void ssmanage(int opt)
{
  ssinfo_t *hdr;

  for(hdr=ssbase_; hdr->next != ssbase_; hdr=hdr->next){
    /* we must not operate on hdr itself, which renders the iterator hdr invalid */
    if(opt == 0) 
      sslistremove_(hdr->next);
    else if(opt == 1)
      ssresize_(&(hdr->next), strlen((char *)(hdr->next+1)), 0);
  }
}

/* copy t to *ps, if ps is not NULL, and return the result 
 * if ps or *ps is NULL, we return a string created from t
 *   *ps is set to the same value if ps is not NULL
 * otherwise, we update the record that corresponds to *ps 
 * */
char *sscpyx(char **ps, const char *t)
{
  ssinfo_t *hdr=NULL;
  size_t size=0;
  char *s, *p;

  if(ps != NULL && (s=*ps) != NULL && (hdr=sslistfind_(s)) == NULL) 
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
    printf("sscpyx: cap=%u, size=%u\n", hdr->size, strlen(*ps));
  }
#endif  
  if((s=ssresize_(&hdr, size, 1)) == NULL) 
    return NULL;
#ifdef SSDBG_
  printf("sscpyx: after resizing, cap=%u\n", hdr->size);
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
  ssinfo_t *hdr;
  size_t size=0;
  char *s, *p;
  const char *q;

  if(ps == NULL || (s=*ps) == NULL || (hdr=sslistfind_(s)) == NULL || t == NULL) 
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
      size, p-s, hdr->size, s, t);
#endif
  if((s=ssresize_(&hdr, size, 1)) == NULL)
    return NULL;
#ifdef SSDBG_
  printf("sscatx: after resize s=%s\naddr: %p\nprev: %p\ncap=%u\n", 
      s, s, *ps, hdr->size);
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

