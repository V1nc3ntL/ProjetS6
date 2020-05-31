/******************************************************************************
			PROJET 	S6
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
*******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <complex.h>
#include <time.h>
#include "tst_func.h"
#include "string.h"
#define MAX 1<<14
#define INIT(a,b) (*a=0,*b=0)

int
main (int argc, char **argv)
{



  struct timespec now, bf;
  FILE *signal, *tf_f, *results;

  clock_gettime (CLOCK_REALTIME, &now);


  results = fopen ("tst_results.txt", "w");


  if (argc == 1)
    no_arg (bf, now, results);
  else
    {
      if (strspn (argv[1], "0123456789") == strlen (argv[1]))
	{
	  tst_one_tf (argv[1], bf, now, results);
	}
    }
  fclose (results);
  return 0;
}
