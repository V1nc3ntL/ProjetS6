/******************************************************************************
                        PROJET  S6
        EL-KHARROUBI    GRAVES  LEFEBVRE        NEZET
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define BUF_SZ 16
#define MAX_SMP_VALUE 1024
FILE *
init_smp (FILE * f, int N)
{
  int i;
  char buf[BUF_SZ];

  for (i = 0; i < N; i++)
    {
      snprintf (buf, BUF_SZ, "%d", i);
      fwrite (buf, sizeof (char), strlen (buf), f);
      fwrite (",", sizeof (char), 1, f);

    }
  return f;
}


 /*!
  * @function   init_smp_rnd
  * @abstract   Initialise un fichier de signal ouvert avec des valeurs comprise entre 0 et MAX_SMP_VALUE
  * @param      f       
  * @result     result
  */
FILE *
init_smp_rnd (FILE * f, int p2, int max)
{
  int i, mx = 1;
  char buf[BUF_SZ];
  mx <<= p2;

  for (i = 0; i < mx; i++)
    {
      snprintf (buf, BUF_SZ, "%d", rand () % max);
      fwrite (buf, sizeof (char), strlen (buf), f);
      fputc ('\n', f);
    }
  return f;
}
