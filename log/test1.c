#include <stdio.h>
#include "log.c"

int main(void)
{
  logfile_t *mylog;
  char msg[1024];

  if ((mylog = log_open("my.log")) == NULL) 
    return 1;
  printf("please write something: ");
  log_printf(mylog, "the input is [%s]\n", 
      fgets(msg, sizeof msg, stdin));
  log_hardflush(mylog);
  printf("log file is hard flushed. please check...");
  getchar();
  log_printf(mylog, "finished.\n");
  log_close(mylog);
  return 0;
}

