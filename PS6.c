/******************************************************************************
FFT²
*******************************************************************************/
#include <stdio.h>
#include <time.h> 
#include <stdlib.h>
#include <math.h>
#define PI              3.1415927
#define ADC_DEF         1024
#define SMP_N       0x4000
#define RE(n,b) *(n+b*2)
#define I(n) *(n+1)
#define IM(n,b) *(n+b*2+1)
void
fftf (short *sig, float *TF, int N, float *twiddles)
{
  int i, j, tmp_N, k, btf, lg2 = 0;
  short *sig_t;
  float tmp_r, tmp_i;
  float *mid, *twid, *base;

  //DFT2

  for (i = N; i; i >>= 1)
    lg2++;

  for (i = 0; i < N; i += 2)
    {
      sig_t = sig + i;
      *(TF + (i << 1)) = *sig_t + *(sig_t + 1);
      *(TF + (i << 1) + 2) = *sig_t - *(sig_t + 1);
    }



  for (tmp_N = 4, i = 1; i < lg2; i++, tmp_N <<= 1)
    {

      for (j = 0; j < N; j += tmp_N)
	{
	  base = TF + (j << 1);
	  mid = base + tmp_N;

	  tmp_N >>= 1;

	  k = (SMP_N) / tmp_N;

	  tmp_r = *(mid);
	  tmp_i = I (mid);

	  I (base) = I (base) + tmp_i;
	  *(mid) = *(base) - tmp_r;
	  I (mid) = I (base) - tmp_i;
	  *(base) = *(base) + tmp_r;

	  for (btf = 1; btf < tmp_N; btf++)
	    {
	      twid = twiddles + k * btf;
	      // Multiplication du twiddle
	      tmp_r = RE (mid, btf) * *(twid) - IM (mid, btf) * I (twid);
	      tmp_i = RE (mid, btf) * I (twid) + IM (mid, btf) ** (twid);

	      RE (mid, btf) = RE (base, btf) - tmp_r;
	      IM (mid, btf) = IM (base, btf) - tmp_i;

	      RE (base, btf) = RE (base, btf) + tmp_r;
	      IM (base, btf) = IM (base, btf) + tmp_i;
	      twid += k;
	    }
	  tmp_N <<= 1;
	}
    }
}

     // Permet d'obtenir les twiddles factors
     // Pour une fft sur b termes
float *
get_twiddles (int b)
{

  int N = b << 1;

  float *twiddles = (float *) malloc (sizeof (float) * b);
  int half = N >> 1, ind;

  for (ind = 0; ind < b; ind += 2)
    {
      *(twiddles + ind) = cos (ind * PI / half);
      *(twiddles + ind + 1) = -sin (ind * PI / half);
    }

  return twiddles;
}

     //Tableau de renversement
static unsigned char br_4[] = {
  0x0, 0x8, 0x4, 0xC, 0x2, 0xA, 0x6, 0xE,
  0x1, 0x9, 0x5, 0xD, 0x3, 0xB, 0x7, 0xF
};

     // Renverse les bits d'un short pour radix 2
int
rv_16_rdx2 (unsigned int j)
{
  return (0x0000ffff | (br_4[(j & 0x000f)] << 12) |
	  (br_4[(j & 0x00f0) >> 4] << 8) |
	  (br_4[(j & 0x0f00) >> 8] << 4) | (br_4[(j & 0x0f000) >> 12]));
}

     // Renverse un tableau de short pour radix 2
void
rvs_16_rdx2 (int *sig)
{
  unsigned int i;
  int tm[SMP_N];
  for (i = 0; i < SMP_N; i++)
    // i renversé
    tm[i] = sig[rv_16_rdx2 (i)];

  for (i = 0; i < SMP_N; i++)
    sig[i] = tm[i];

}

int
main ()
{
  unsigned int i, N = SMP_N, b = N << 1;

  /* Tableau de complexes de taille N */
  // float *TF = (float *) calloc (b, sizeof (float) * b);
  float *TF = (float *) calloc (b, sizeof (float) * b);
  short *sig = malloc (sizeof (short) * (N));
  float *twiddles = get_twiddles (N);
	struct timespec now,bf;
  if (!sig || !twiddles || !TF)
    return 0;

  printf ("TF sur %d termes\n", N);
  printf ("\n\tBesoins mémoires :");
  printf ("\nTF :\t %d o \n", b * 4);
  printf ("SIG : \t%d o \n", N * 2);

  // Tablieau renversant le tableau du signal

  for (i = 0; i < N; i++)
    *(sig + i) = i;

//  rvs_16_rdx2 (sig);
    timespec_get( &bf, TIME_UTC );
 	fftf (sig, TF, N, twiddles);

    timespec_get( &now, TIME_UTC );
    printf("\nTemps écoulé : %lf ms\n" ,(now.tv_nsec-bf.tv_nsec)/1e3);

  for (int i = 0; i < 8; i += 2)
    {
      printf ("R %-2.2f\t\t", *(TF + i));
      printf ("%-2.2f I\n", I (TF + i));
    }

  free (TF);
  free (twiddles);
  free (sig);
  return 0;
}
