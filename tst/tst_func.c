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
#include <string.h>
#include "br_tst.h"
#include "tst_func.h"
#include "fft_tst.h"
#include "fft_cmp_tst.h"
#include "string.h"
#include "smp_tst.h"
#include "cmp_tf.h"
#include "../src/fft.h"
#define BUF_SZ 64
#define MAX_POW 14
#define MAX 1<<MAX_POW
#define MAX_SMP_VAL 1024
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
  int tmp_N;
  char *mess[2] = { "TF sur ",
    " termes\n\t\t\tSignal\tTF\tTwiddles\t\nBesoins mémoires :\t"
  };  
  fputs (*mess, res);
  if (N > (1<<13) ){
	tmp_N = (N>>10);
  }else{
      tmp_N = N;
  }


  snprintf (buf, BUF_SZ, "%d", N); 
  fwrite (buf, sizeof (char), strlen (buf), res);
  fputs (*(mess + 1), res);
  snprintf (buf, BUF_SZ, "%d", (int) (tmp_N * sizeof (short)));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fputc ('\t', res);

  snprintf (buf, BUF_SZ, "%d", (int) (tmp_N * sizeof (float) * 2));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fwrite ("\t\t", sizeof (char), 2, res);
  snprintf (buf, BUF_SZ, "%d", (int) (tmp_N * sizeof (float) * 2));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fputc ('\t', res);
  if (N > (1<<13) ){
  fputc('k', res);
  }

  fputc('o', res);
}

void
write_a_tf_cmp (int N, FILE * res)
{
  char buf[BUF_SZ];
  char *mess[2] = { "TF sur ",
    " termes\n\t\t\tSignal\tTF\tTF complexe\tTwiddles\t\nBesoins mémoires :\t"
  };
  fputs (*mess, res);
  snprintf (buf, BUF_SZ, "%d", N);
  fwrite (buf, sizeof (char), strlen (buf), res);
  fputs (*(mess + 1), res);
  snprintf (buf, BUF_SZ, "%d", (int) (N * sizeof (short)));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fputc ('\t', res);

  snprintf (buf, BUF_SZ, "%d", (int) (N * sizeof (float) * 2));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fputc ('\t', res);
  snprintf (buf, BUF_SZ, "%d", (int) (N * sizeof (complex)));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fwrite ("\t\t", sizeof (char), 2, res);
  snprintf (buf, BUF_SZ, "%d", (int) (N * sizeof (float) * 2));
  fwrite (buf, sizeof (char), strlen (buf), res);
  fputs ("\to", res);
  fputc ('\t', res);
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
    { "\nTemps execution :\n", "\tRenversement 2 :\t", "\tRenversement 4 :\t",
    "\tRADIX 2 :\t\t", "\tRADIX 4 :\t\t", "\tFFTW :\t\t\t"
  };

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

      fputs (*(mess + 1), res);
      print_time (res, time);

      fputs (*(mess + 3), res);
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

	  fputs (*(mess + 2), res);
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
  int N, i;
  /* Tableau de complexes de taille N */
  float *TF;
  float time;
  /* Tableau du signal */
  short *sig;
  /* Tableau des twiddles */
  float *twiddles;

  char state = 0;

  unsigned int *cmp_mul = (unsigned int *) calloc (1, sizeof (unsigned int)),
    *cmp_add = (unsigned int *) calloc (1, sizeof (unsigned int));
  char *mess[6] =
    { "Temps execution :\n", "\tRenversement 2 :\t", "\tRenversement 4 :\t",
    "\tRADIX 2 :\t\t", "\tRADIX 4 :\t\t", "\tFFTW :\t\t\t"
  };

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

  fputs (*(mess + 1), res);
  print_time (res, time);

  fputs (*(mess + 3), res);
  time = fftf_tst (fftf_rdx2_tst, N, sig, TF, twiddles, &bf, &now, cmp_add,
		   cmp_mul);
  print_time (res, time);
  if (!(atoi (number) & 0x1))
    {
      *cmp_mul = 0;
      *cmp_add = 0;

      for (i = 0; i < N; i++)
	*(sig + i) = i;

      fputs (*(mess + 2), res);
      time = br_tst (rvs_16_rdx4, N, sig, &bf, &now);
      print_time (res, time);
      fputs (*(mess + 4), res);

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

// 


  free (TF);
  free (twiddles);
  free (sig);
  free (cmp_add);
  free (cmp_mul);
  return 1;
}


/*!
 * @function	init
 * @abstract	Créé une fonction pour obtenir un fichier de signal formaté
 * @param	argv2	le fichier de signal
 * @result	result
 */
int
init (char *argv2)
{
  int tmp;
  char buf_c[BUF_SZ];
  FILE *signal;
  strcpy (buf_c, argv2);
  if (strspn (buf_c, "0123456789") == strlen (buf_c))
    {
      tmp = (int) strtol (buf_c, NULL, 10);
      if (tmp > MAX_POW)
	{
	  printf ("Taille maximum dépassée\n");
	  return 0;
	}
      else
	{
	  signal = fopen ("signal.txt", "w");
	  init_smp_rnd (signal, tmp, MAX_SMP_VAL);
	}
      fclose (signal);
      return 1;
    }
  printf ("Saisissez une puissance de 2 en 2nd argument\n");
  return 0;
}

/*!
 * @function    ml_compare
 * @abstract    Compare les algorithmes avec les résultats de matlab stockés dans un fichier
 * @param       sig_name chaîne du nom du fichier de signal     
 * @param       ml_name_name chaîne du nom du fichier de fft matlab     
 * @result      result
 */
int
ml_compare (char *sig_name, char *ml_name,FILE* res,struct timespec bf, struct timespec now)
{
  int N, i = 0, tmp;
  char buf_c[BUF_SZ];
  FILE *signal, *tf_f;
  short *sig;
  float time,tmp_f, *tf, *tf_ml, *twi;

  char *mess[6] =
    { "\nTemps execution :\n", "\tRenversement 2 :\t", "\tRenversement 4 :\t",
    "\tRADIX 2 :\t\t", "\tRADIX 4 :\t\t"
  };
  signal = fopen (sig_name, "r");
  strcpy (buf_c, ml_name);
  tf_f = fopen (buf_c, "r");

  while ((tmp = fgetc (signal)) != EOF)
    {
      if (tmp == '\n')
	i++;
    }
  rewind (signal);

  if (i & 1)
    {
      fclose (signal);
      printf ("Le nombre de terme n'est pas pair ! ");
      return 0;
    }

  N = i;

   write_a_tf ( N,res);


  tf = (float *) calloc (N << 1, sizeof (float) * (N << 1));
  tf_ml = (float *) calloc (N << 1, sizeof (float) * (N << 1));
  twi = get_twiddles (N);
  sig = (short *) malloc (sizeof (short) * (N));

  if (!tf_ml||!tf||!twi || !sig || !tf)
    return 0;


  for (i = 0; i < (N); i++)
    fscanf (tf_f, "%f", tf_ml + i * 2);

  for (i = 0; i < N; i++)
    fscanf (tf_f, "%f", tf_ml + i * 2 + 1);


  fclose (tf_f);
  for (i = 0; i < N; i++)
    fscanf (signal, "%hd", sig + i);


  fputs (*(mess ), res);
  fputs (*(mess + 1), res);
  time = br_tst (rvs_16_rdx2,  N, sig, &bf,&now);
      print_time(res,time);

  time = fftf_tst (fftf_rdx2_tst,  N, sig,tf,twi,&bf,&now,0,0);
  fputs (*(mess + 3), res);
      print_time(res,time);


  tmp_f = cmp_tf (tf, tf_ml, N);

  if (tmp_f != 0)
    fprintf (res,"\tErreur Radix 2:\t\t%e\n\n", tmp_f);

  rewind (signal);
  for (i = 0; i < N; i++)
    fscanf (signal, "%hd", sig + i);
  tmp = N;

  i = 0;
  while (tmp >>= 1)
    i++;

  if (!(i & 0x1))
    {
      fputs (*(mess + 2), res);
      time = br_tst (rvs_16_rdx4,  N, sig, &bf,&now);
      print_time(res,time);
      free (tf);
      tf = (float *) calloc (N << 1, sizeof (float) * (N << 1));
      if(!tf)
	  return 0;

      time = fftf_tst (fftf_rdx4_tst,  N, sig,tf,twi,&bf,&now,0,0);
      fputs (*(mess + 4), res);
      print_time(res,time);

      tmp_f = cmp_tf (tf, tf_ml, N);

      if (tmp_f != 0)
	fprintf (res,"\tErreur Radix 4:\t\t%e\n\n", tmp_f);

    }
  fclose (signal);

  free (tf);
  free (twi);
  free (sig);

  return 1;
}
