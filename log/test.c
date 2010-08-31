#include <stdio.h>
#include "log.h"
#include "ss.h"

int main(void)
{
  logfile_t *mylog;
  char *msg = NULL;

  if ((mylog=logopen("my.log")) == NULL) 
    return 1;
  printf("please write something: ");
  logprintf(mylog, "the input is [%s]\n", 
      ssfgets(msg, NULL, stdin));
  loghardflush(mylog);
  printf("log file is hard flushed. please check...");
  getchar();
  logprintf(mylog, "finished.\n");
  logclose(mylog);
  return 0;
}

