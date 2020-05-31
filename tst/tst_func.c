/******************************************************************************
			PROJET 	S6
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
*******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <complex.h>
#include <time.h>
#include <fftw3.h>
#include "br_tst.h"
#include "tst_func.h"
#include "fft_tst.h"
#include "fft_cmp_tst.h"
#include "string.h"
#define MAX 1<<14
#define BUF_SZ 64
#define INIT(a,b) (*a=0,*b=0)
float
br_tst (void (*br) (short *, int), int N, short *sig, struct timespec *bf,
	struct timespec *now)
{
  timespec_get (bf, TIME_UTC);
  br (sig, N);
  timespec_get (now, TIME_UTC);
  return (now->tv_nsec - bf->tv_nsec) / 1e3;
}


float
fftf_tst (void (*fft)
	  (short *, float *, int, float *, unsigned int *, unsigned int *),
	  int N, short *sig, float *TF, float *twiddles, struct timespec *bf,
	  struct timespec *now, unsigned int *add, unsigned int *mul)
{
  timespec_get (bf, TIME_UTC);
  fft (sig, TF, N, twiddles, add, mul);
  timespec_get (now, TIME_UTC);
  return (now->tv_nsec - bf->tv_nsec) / 1e3;
}

float
fftf_cmp_tst (void (*fft)
	      (short *, complex *, int, complex *, unsigned int *,
	       unsigned int *), int N, short *sig, complex * TF,
	      complex * twiddles, struct timespec *bf, struct timespec *now,
	      unsigned int *add, unsigned int *mul)
{
  timespec_get (bf, TIME_UTC);
  fft (sig, TF, N, twiddles, add, mul);
  timespec_get (now, TIME_UTC);
  return (now->tv_nsec - bf->tv_nsec) / 1e3;
}

float
fftw_tst (int N, short *sig, float *TF, struct timespec *bf,
	  struct timespec *now)
{
  fftw_complex *in, *out;
  int i;
  fftw_plan *p = (fftw_plan *) fftw_malloc (sizeof (fftw_plan));
  in = fftw_alloc_complex (N);
  out = fftw_alloc_complex (N);
  for (i = 0; i < N; i++)
    in[i] = sig[i];

  *p = fftw_plan_dft_1d (N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  timespec_get (bf, TIME_UTC);
  fftw_execute (*p);
  timespec_get (now, TIME_UTC);
  for (i = 0; i < N; i += 2)
    {
      TF[i] = creal (out[i]);
      TF[i + 1] = cimag (out[i]);
    }
  fftw_free (in);
  fftw_free (out);
  fftw_free (p);
  return (now->tv_nsec - bf->tv_nsec) / 1e3;
}

void
write_a_tf (int N, FILE * res)
{
  char buf[BUF_SZ];
  char *mess[2] =
    { "TF sur ",
" termes\n\t\t\tSignal\tTF\tTF complexe\tTwiddles\t\nBesoins mémoires :\t" };
  fwrite (*mess, sizeof (char), sizeof (*mess) - 1, res);
  snprintf (buf, BUF_SZ, "%d", N);
  fwrite (buf, sizeof (char), strlen (buf), res);
  fwrite (*(mess + 1), sizeof (char), strlen (*(mess + 1)), res);
  snprintf (buf, BUF_SZ, "%d", (int) (N * sizeof (short)));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fwrite ("\t", sizeof (char), 1, res);
  snprintf (buf, BUF_SZ, "%d", (int) (N * sizeof (float) * 2));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fwrite ("\t", sizeof (char), 1, res);
  snprintf (buf, BUF_SZ, "%d", (int) (N * sizeof (complex)));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fwrite ("\t\t", sizeof (char), 2, res);
  snprintf (buf, BUF_SZ, "%d", (int) (N * sizeof (float) * 2));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fwrite ("\to", sizeof (char), 2, res);
  fwrite ("\n", sizeof (char), 1, res);
}

void
print_time (FILE * res, float time)
{
  char buf[BUF_SZ];
  snprintf (buf, BUF_SZ, "%+.3e", time);
  fwrite (buf, sizeof (char), strlen (buf), res);
  fwrite ("\t ms\n", sizeof (char), 5, res);
}

int
no_arg (struct timespec bf, struct timespec now, FILE * res)
{
  char state = 0;
  int N, i;
  /* Tableau de complexes de taille N */
  float *TF;
  float time;
  complex *TF_c;
  /* Tableau du signal */
  short *sig;
  unsigned int *cmp_mul = (unsigned int *) calloc (1, sizeof (unsigned int)),
    *cmp_add = (unsigned int *) calloc (1, sizeof (unsigned int));
  /* Tableau des twiddles */
  float *twiddles;
  complex *twiddles_c;
  char *mess[6] =
    { "Temps execution :\n", "\tRenversement 2 :\t", "\tRenversement 4 :\t",
"\tRADIX 2 :\t\t", "\tRADIX 4 :\t\t", "\tFFTW :\t\t\t" };

  for (N = 16; N <= MAX; N <<= 1)
    {
      twiddles = get_twiddles (N);
      twiddles_c = get_twiddles_cmp (N);
      sig = calloc (N, sizeof (short) * (N));
      TF = (float *) calloc (N << 1, sizeof (float) * (N << 1));
      TF_c = (complex *) calloc (N, sizeof (fftw_complex) * (N));
      if (!sig || !twiddles || !TF || !cmp_add || !cmp_mul || !TF_c
	  || !twiddles_c)
	{
	  return 0;
	}
      write_a_tf (N, res);

      for (i = 0; i < N; i++)
	*(sig + i) = i;

      fwrite (*mess, sizeof (char), strlen (*mess), res);
      time = br_tst (rvs_16_rdx2, N, sig, &bf, &now);

      fwrite (*(mess + 1), sizeof (char), strlen (*(mess + 1)), res);
      print_time (res, time);

      fwrite (*(mess + 3), sizeof (char), strlen (*(mess + 3)), res);
      time =
	fftf_tst (fftf_rdx2_tst, N, sig, TF, twiddles, &bf, &now, cmp_add,
		  cmp_mul);
      print_time (res, time);
      if (!state)
	{
	  *cmp_mul = 0;
	  *cmp_add = 0;
	  for (i = 0; i < N; i++)
	    *(sig + i) = i;

	  fwrite (*(mess + 2), sizeof (char), strlen (*(mess + 2)), res);
	  time = br_tst (rvs_16_rdx4, N, sig, &bf, &now);
	  print_time (res, time);
	  fwrite (*(mess + 4), sizeof (char), strlen (*(mess + 4)), res);

	  time = fftf_tst (fftf_rdx4_tst, N, sig, TF, twiddles, &bf, &now,
			   cmp_add, cmp_mul);
	  print_time (res, time);
	}
      *cmp_mul = 0;
      *cmp_add = 0;

      for (i = 0; i < N; i++)
	*(sig + i) = i;

      fwrite (*(mess + 5), sizeof (char), strlen (*(mess + 5)), res);
      time = fftw_tst (N, sig, TF, &bf, &now);

      print_time (res, time);

      state = !state;
      free (TF);
      free (TF_c);
      free (twiddles);
      free (twiddles_c);
      free (sig);
    }
  free (cmp_add);
  free (cmp_mul);
  return 1;
}

int
tst_one_tf (char *number, struct timespec bf, struct timespec now, FILE * res)
{
  char state = 0;
  int N, i;
  /* Tableau de complexes de taille N */
  float *TF;
  float time;
  /* Tableau du signal */
  short *sig;
  unsigned int *cmp_mul = (unsigned int *) calloc (1, sizeof (unsigned int)),
    *cmp_add = (unsigned int *) calloc (1, sizeof (unsigned int));
  /* Tableau des twiddles */
  float *twiddles;
  char *mess[6] =
    { "Temps execution :\n", "\tRenversement 2 :\t", "\tRenversement 4 :\t",
"\tRADIX 2 :\t\t", "\tRADIX 4 :\t\t", "\tFFTW :\t\t\t" };

  N = (int) strtol (number, NULL, 10);
  N = 1 << (N);
  printf ("%d points\n", N);

  twiddles = get_twiddles (N);
  sig = calloc (N, sizeof (short) * (N));
  TF = (float *) calloc (N << 1, sizeof (float) * (N << 1));
  if (!sig || !twiddles || !TF || !cmp_add || !cmp_mul)
    {
      return 0;
    }
  write_a_tf (N, res);

  for (i = 0; i < N; i++)
    *(sig + i) = i;

  fwrite (*mess, sizeof (char), strlen (*mess), res);
  time = br_tst (rvs_16_rdx2, N, sig, &bf, &now);

  fwrite (*(mess + 1), sizeof (char), strlen (*(mess + 1)), res);
  print_time (res, time);

  fwrite (*(mess + 3), sizeof (char), strlen (*(mess + 3)), res);
  time = fftf_tst (fftf_rdx2_tst, N, sig, TF, twiddles, &bf, &now, cmp_add,
		   cmp_mul);
  print_time (res, time);
  if (!(atoi (number) & 0x1))
    {
      *cmp_mul = 0;
      *cmp_add = 0;

      for (i = 0; i < N; i++)
	*(sig + i) = i;

      fwrite (*(mess + 2), sizeof (char), strlen (*(mess + 2)), res);
      time = br_tst (rvs_16_rdx4, N, sig, &bf, &now);
      print_time (res, time);
      fwrite (*(mess + 4), sizeof (char), strlen (*(mess + 4)), res);

      time = fftf_tst (fftf_rdx4_tst, N, sig, TF, twiddles, &bf, &now,
		       cmp_add, cmp_mul);
      print_time (res, time);
    }
  *cmp_mul = 0;
  *cmp_add = 0;

  for (i = 0; i < N; i++)
    *(sig + i) = i;

  fwrite (*(mess + 5), sizeof (char), strlen (*(mess + 5)), res);
  time = fftw_tst (N, sig, TF, &bf, &now);

  print_time (res, time);

  state = !state;
  free (TF);
  free (twiddles);
  free (sig);
  free (cmp_add);
  free (cmp_mul);
  return 1;
}
