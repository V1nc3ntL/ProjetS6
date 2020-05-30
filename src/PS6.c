/******************************************************************************
			PROJET 	S6
			FFT
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
FFT : radix 2 et 4
*******************************************************************************/


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "bit_reversal.h"
#include "fft.h"

#define PI 3.1415927
#define ADC_DEF 1024


#define SMP_N 1<<14 

int
main (void)
{


  int i, N = SMP_N, b = N << 1;

/* Tableau de complexes de taille N */
  float *TF;

  short *sig = calloc (N, sizeof (short) * (N));
  TF = (float *) calloc (b, sizeof (float) * b);
  float *twiddles = get_twiddles (N);



  struct timespec now, bf;


  if (!sig || !twiddles || !TF)
    {
      printf ("Pas assez de mémoire\n");
      return 0;
    }


  printf ("TF sur %d termes\n", N);


  printf ("\n\tBesoins mémoires :");


  printf ("\nTF :\t %d o \n", b * 4);


  printf ("SIG : \t%d o \n", N * 2);



  for (i = 0; i < N; i++)
    *(sig + i) = i;



  rvs_16_rdx2 (sig,N);
/*
  for (i = 0; i < N; i++)
    {
      if (!(i % 4) && i)
	printf ("\n");
      printf ("%d,", sig[i]);
    }
*/
  timespec_get (&bf, TIME_UTC);
 fftf_rdx2 (sig, TF, N, twiddles);
  timespec_get (&now, TIME_UTC);
  printf ("\nRADIX 2\nTemps écoulé :\t %lf ms\n",
	  (now.tv_nsec - bf.tv_nsec) / 1e3);

  for (i = 0; i < N; i++)
    *(sig + i) = i;

  timespec_get (&bf, TIME_UTC);

  rvs_16_rdx4 (sig, N);

  fftf_rdx4 (sig, TF, N, twiddles);

  timespec_get (&now, TIME_UTC);

  printf ("\n");
  printf ("\nRADIX 4\nTemps écoulé :\t %*lf ms\n", 4,
	  (now.tv_nsec - bf.tv_nsec) / 1e3);
  /*
for(i = 0 ; i < N ; i++)	{
      printf ("%d R %+-2.2f\t\t", i, *(TF+i*2));
      printf ("%+-2.2f I\n", *(TF + i*2+ 1));
}*/
  free (TF);
  free (twiddles);
  free (sig);
  return 0;
}
