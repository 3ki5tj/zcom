#include <stdio.h>
#include "log.c"
#include "ss.c"

int main(void)
{
  logfile_t *mylog;
  char *msg = NULL;

  if ((mylog = log_open("my.log")) == NULL) 
    return 1;
  printf("please write something: ");
  log_printf(mylog, "the input is [%s]\n", 
      ssfgets(msg, NULL, stdin));
  log_hardflush(mylog);
  printf("log file is hard flushed. please check...");
  getchar();
  log_printf(mylog, "finished.\n");
  log_close(mylog);
  return 0;
}

