#ifndef SAFESTRING
#define SAFESTRING 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ss.h"

#define SSMINSIZ (1u<<8)

typedef struct tag_ssinfo{
  size_t size;
  struct tag_ssinfo *prev, *next;
  size_t mark;
}ssinfo_t;

#define SSMARK (size_t)0x1f3d5b79 /* signature mark */

static struct tag_ssinfo ssbase_[1]={{0u,NULL,NULL,SSMARK}};

#define sserror_ printf

/* calculate memory required to accommendate a size-byte string */
static size_t sscalcsize_(size_t size, int hdr)
{
  size = (size/SSMINSIZ+1) * SSMINSIZ;
  if(hdr) size+=sizeof(ssinfo_t);
  return size;
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

/* (re)allocate memory for *phdr, update list, return the new string */
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
        *phdr = hdr;
      }
      hdr->size = size;
      hdr->mark = SSMARK;
    }
  }
  return (char *)(hdr+1);
}

static void ssdump_(ssinfo_t *hdr)
{
  printf("header: %p, prev: %p, next: %p, base: %p, size: %6d, string: %p\n", 
      hdr, hdr->prev, hdr->next, ssbase_, hdr->size, (char *)hdr+sizeof(*hdr));
}

/* return a newly created string with its content copied from t */
char *ssnew(const char *t)
{
  ssinfo_t *hdr=NULL;
  size_t size;
  char *s;

  size = t ? strlen(t) : 0;
  if((s=ssresize_(&hdr, size, 1)) == NULL)
    return NULL;
  else if(t)
    return strcpy(s, t);
  else{
    s[0]='\0';
    return s;
  }
}

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

char *sscopy_x(char **ps, const char *t)
{
  ssinfo_t *hdr;

  if(ps == NULL || t == NULL)
    return NULL;

  if(*ps == NULL) return *ps = ssnew(t);

  if((hdr=sslistfind_(*ps)) == NULL) 
    return NULL;

  if((*ps=ssresize_(&hdr, strlen(t)+1, 1)) == NULL) 
    return NULL;
  else
    return strcpy(*ps, t);
}

char *sscat_x(char **ps, const char *t)
{
  ssinfo_t *hdr;
  size_t size;

  if(ps == NULL || *ps == NULL || (hdr=sslistfind_(*ps)) == NULL || t == NULL) 
    return NULL;

  size = strlen(*ps);
  size += strlen(t)+1;
  if((*ps=ssresize_(&hdr, size, 1)) == NULL)
    return NULL;
  else 
    return strcat(*ps, t);
}

#endif

