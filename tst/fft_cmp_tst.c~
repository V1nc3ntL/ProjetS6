#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define PI 3.1415927
#define X3(n) ((n<<1) + n)




void
fftf_rdx2_cmp_tst (short *sig, complex * TF, int N, complex * twiddles,
		   unsigned int *cmp_add, unsigned int *cmp_mul)
{
  int i, j, tmp_N, k, btf, lg2 = 0;
  short *sig_t;
  complex tmp;
  complex *mid, *twid, *base;
//DFT2
  for (i = N; i; i >>= 1)
    lg2++;

// Toutes les DFT 2 sont effectuées avec des additions réelles
  for (i = 0; i < N; i += 2)
    {
      sig_t = sig + i;
      *(TF + i) = *sig_t + *(sig_t + 1);
      *(TF + i + 1) = *sig_t - *(sig_t + 1);
    }


  for (tmp_N = 4, i = 1; i < lg2; i++, tmp_N <<= 1)
    {
      for (j = 0; j < N; j += tmp_N)
	{
	  /*Adresse de X[0+j] 
	     et de X[N/2+j]
	   */

	  k = (N) / tmp_N;
	  tmp_N >>= 1;
	  base = TF;
	  mid = base + tmp_N;
	  tmp = *(mid);

	  *mid = *base - tmp;
	  *base = *base + tmp;

	  *cmp_add += 2;

	  for (btf = 1; btf < tmp_N; btf++)
	    {
	      //W btf 
	      twid = twiddles + k * btf;
	      *cmp_mul += 1;

	      tmp = *(mid + btf) * *(twid);

	      *(mid + btf) = *(base + btf) - tmp;
	      *(base + btf) = *(base + btf) + tmp;
	      *cmp_add += 2;

	      twid += k;
	    }
	  tmp_N <<= 1;
	}
    }
}

void
fftf_rdx4_cmp_tst (short *sig, complex * TF, int N, complex * twiddles,
		   unsigned int *cmp_add, unsigned int *cmp_mul)
{


  int lvl, i, j, k, tmp_N, btf, lg4 = 0;

  complex bf[4];

  complex *twid, *pos[4];

  for (i = N; i; i >>= 2)
    lg4++;
  lg4--;

  //W0
  *TF = (complex) * (sig) + *(sig + 1) + *(sig + 2) + *(sig + 3);
  //W1
  *(TF + 1) =
    (complex) * (sig) - I * *(sig + 1) - *(sig + 2) + I ** (sig + 3);
  //W2
  *(TF + 2) = (complex) * sig - *(sig + 1) + *(sig + 2) - *(sig + 3);
  //W3
  *(TF + 3) = (complex) * sig + I ** (sig + 1) - *(sig + 2) - I ** (sig + 3);
//FFT 4
  for (i = 4; i < N; i += 4)
    {
      *pos = TF + i;
      //W0
      **pos = (complex) * (sig) + *(sig + 1) + *(sig + 2) + *(sig + 3);
      //W1
      *pos[1] =
	(complex) * (sig) - I * *(sig + 1) - *(sig + 2) + I ** (sig + 3);
      //W2
      *pos[2] = (complex) * sig - *(sig + 1) + *(sig + 2) - *(sig + 3);
      //W3
      *pos[3] =
	(complex) * sig + I ** (sig + 1) - *(sig + 2) - I ** (sig + 3);
    }


//FFT(4^(lvl+1)) 
  for (tmp_N = 16, lvl = 1; lvl < lg4; lvl++, tmp_N <<= 2)
    {
      k = N / tmp_N;

//Sur l'ensemble du signal
      for (j = 0; j < N; j += tmp_N)
	{
	  //Position de FFT courante
	  *pos = TF + j;
	  //N/4 (TF*2 termes)
	  tmp_N >>= 2;

	  //Récupération N/4,N/2,3N/4
	  for (i = 1; i < 4; i++)
	    pos[i] = pos[i - 1] + (tmp_N);

	  //Calcul du premier papillon
	  for (i = 0; i < 4; i++)
	    bf[i] = *pos[i];

	  //X(0) 
	  **pos = (complex) * (bf) + *(bf + 1) + *(bf + 2) + *(bf + 3);

	  //X(N/4)
	  *pos[1] =
	    (complex) * (sig) - I * *(sig + 1) - *(sig + 2) + I ** (sig + 3);
	  //X(N/2)ne peut être obtenu par symétrie
	  *pos[2] = (complex) * sig - *(sig + 1) + *(sig + 2) - *(sig + 3);
	  *cmp_add += 12;

	  for (btf = 1; btf < tmp_N; btf++)
	    {

	      //Haut du papillon
	      *pos = TF + j + btf;
	      *bf = *pos[0];
	      //Récupération N/4,N/2,3N/4
	      for (i = 1; i < 4; i++)
		pos[i] = pos[i - 1] + (tmp_N);

	      for (i = 1; i < 4; i++)
		{
		  //W (i * btf) * N i/4
		  twid = twiddles + (i * (btf) * k * 2);
		  *(bf + i) = *pos[i] ** twid;
		  *cmp_mul += 1;
		}
	      //X(4k) 
	      **pos = *(bf) + *(bf + 1) + *(bf + 2) + *(bf + 3);
	      //X(4k+1) 
	      *pos[1] =
		*bf + I * (*(bf + 1)) - *(bf + 2) - I * cimag (*(bf + 3));
	      *cmp_add += 8;
	    }
	  //Les DFTS réelles sont symétriques
	  *pos = TF + j + (tmp_N << 1);
	  for (i = 1; i < tmp_N; i++)
	    {
	      *pos[i] = conj (*(*(pos - i)));
	    }

	  // N courant
	  tmp_N <<= 2;
	}
    }
}

// Pour une fft sur N termes
complex *
get_twiddles_cmp (int N)
{
  int i;
  complex *twiddles = (complex *) malloc (sizeof (complex) * N);

  for (i = 0; i < N; i++)
    {
      //e -i2pi/N 
      *(twiddles + i) = cexp (-i * 2 * PI / N);
    }
  return twiddles;
}
