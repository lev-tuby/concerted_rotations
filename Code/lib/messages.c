#include <stdlib.h>
#include <stdio.h>
#include "messages.h"


void failed (char message[])
{
  fprintf (stderr,"Error. \n FAILED *** %s ***\n \n", message);
  abort ();
  //exit (1);
}

void checkpoint(char message[]){

printf ("\n CHECKPOINT *** %s ***\n \n", message);

}

void error(char message[], char *path, int line)
{
    fprintf(stderr,"ERROR:\n%s\nIn file:%s line:%d\n", message, path, line);
    fflush(stderr);
    abort();
}

void warning(char message[], char *path, int line)
{
    if (global_warn < global_max_warn)
    {
        fprintf(stderr, "WARNING:\n%s\n", message);
        fflush(stderr);
        global_warn++;
    }
    else
    {
        fprintf(stderr, "WARNING:\n%s\nIn file:%s line:%d\n", message, path, line);
        fflush(stderr);
        abort();
    }
}

void info(char message[], char *path, int line)
{
    fprintf(stdout, "INFO:\n%s\n", message);
    fflush(stdout);
}

