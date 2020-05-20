/******************************************************************************
FFT²
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI              3.1415927
#define ADC_DEF         1024
#define SMP_N       1<<16
#define RE(n,b) *(n+b*2)
#define I(n) *(n+1)
#define IM(n,b) *(n+b*2+1)
void
fftf (int *sig, float *TF, int N, float *twiddles)
{
  int i, j, tmp_N, k, btf, lg2 = 0;
  int * sig_t;
  float tmp_r, tmp_i;
  float *mid, *twid, *base;

  //DFT2

  for (i = 0; i < N; i += 2)
    {
      sig_t = sig + i;
      *(TF + (i << 1)) = *sig_t + *(sig_t + 1);
      *(TF + (i << 1) + 2) = *sig_t - *(sig_t + 1);
    }

  for (i = N; i >>= 1; i)
    lg2++;

  printf ("LOG 2 : %d\n", lg2);

  for (tmp_N = 4, i = 1; i < lg2; i++, tmp_N <<= 1)
    {

      for (j = 0; j < N; j += tmp_N)
	{
	  base = TF + (j << 1);
	  mid = base + tmp_N;

	  tmp_N>>=1;

	  k = (SMP_N) / tmp_N ;
	  
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
		printf("\nTW : %f",*twid);
		printf("\nTWIM %f\n",*(twid+1));
	      tmp_r = RE (mid, btf) * *(twid) - IM (mid, btf) * I (twid);
	      tmp_i = RE (mid, btf) * I (twid) + IM (mid, btf) ** (twid);

	      RE (mid, btf) = RE (base, btf) - tmp_r;
	      IM (mid, btf) = IM (base, btf) - tmp_i;

	      RE (base, btf) = RE (base, btf) + tmp_r;
	      IM (base, btf) = IM (base, btf) + tmp_i;
	      twid += k;

	    }
	    tmp_N<<=1;
	 
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
  return ((br_4[(j & 0x000f)] << 12) |
	  (br_4[(j & 0x00f0) >> 4] << 8) |
	  (br_4[(j & 0x0f00) >> 8] << 4) | (br_4[(j & 0x0f000) >> 12]));
}

     // Renverse un tableau de short pour radix 2
void
rvs_16_rdx2 (int *sig)
{
  int i;
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
    int i, N = 1<<14 , b = N << 1;
  int *sig = malloc (sizeof (int) * (N));
  //int * sig=(int*) calloc(SMP_N,sizeof(int)*SMP_N) ;
 
  /* Tableau de complexes de taille N*2 */
  float *TF = (float *) calloc (b, sizeof (float) * b);
  float *twiddles = get_twiddles (SMP_N);
  // Tableau renversant le tableau du signal
  //  rvs_16_rdx2 (sig);

 	for(i=0; i < N; i++){
		*(sig+i) = i ;
	}

  fftf (sig, TF, N, twiddles);
  for (int i = 0; i < b; i += 2)
    {
      printf ("R %-2.2f\t\t", *(TF + i));
      printf ("%-2.2f I\n", I (TF + i));
    }

  free (TF);
  free (twiddles);
  free (sig);
  return 0;
}
